"""Fast engine: everything in C (gaiascan.so) — libdeflate + scan + write.

The C kernel decompresses each gzip file with libdeflate and scans it in a
single pass (extracting source_id and the min/max of the positive BP/RP fluxes)
across all files in parallel with the GIL released, then computes the
percentage change, filters > 100%, sorts, and writes the challenge CSV — all in
C. Python only marshals the file paths and the output path, so there is no
Polars/NumPy import cost and no per-row Python work.

Robustness: the .so is compiled at image build; if it's missing at run time
(e.g. a dev bind-mount shadows it) it is compiled on first use, and if that
fails the module raises so RunScript can fall back to the Polars engine.
"""
import ctypes
import glob
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
# The kernel is one fused decompress+scan pass per file. The job is bound by the
# largest (unsplittable) gzip stream, so once we schedule that file first, extra
# threads only add memory-bandwidth contention on the parallel decompress —
# ~0.75x cores measured fastest. Passed into the kernel (omp_set_num_threads)
# rather than via OMP_NUM_THREADS, since libgomp is already initialized by IRIS.
# Use the affinity mask (honors Docker --cpuset-cpus / cgroup limits) rather than
# os.cpu_count(), which reports the whole host and would oversubscribe the kernel
# 6:1 when the container is pinned to a subset of cores (e.g. 4 of 32).
try:
    _NCPU = len(os.sched_getaffinity(0))
except AttributeError:  # non-Linux
    _NCPU = os.cpu_count() or 8
NTHREADS = max(4, int(_NCPU * 3 // 4))
SO = os.path.join(HERE, "gaiascan.so")
SRC = os.path.join(HERE, "gaiascan.c")
IN_DIR = "/home/irisowner/dev/data/in"
OUT = "/home/irisowner/dev/data/out/challenge_output.csv"
SLOT = 8000  # max rows reserved per file (files hold ~5000)


def _find_libdeflate():
    for p in ("/usr/lib/aarch64-linux-gnu/libdeflate.so.0",
              "/usr/lib/x86_64-linux-gnu/libdeflate.so.0",
              "/usr/lib/libdeflate.so.0"):
        if os.path.exists(p):
            return p
    raise RuntimeError("libdeflate.so.0 not found")


def _ensure_so():
    if os.path.exists(SO):
        return
    # -march=native is safe here: we compile on the same machine that runs it.
    # Fall back to a portable build if a native build fails.
    base = ["gcc", "-fopenmp", "-fPIC", "-shared", SRC, _find_libdeflate(), "-o", SO]
    native = ["-O3", "-march=native", "-funroll-loops"]
    if subprocess.run(base[:1] + native + base[1:]).returncode != 0:
        subprocess.run(base[:1] + ["-O3"] + base[1:], check=True)


_ensure_so()
_LIB = ctypes.CDLL(SO)
_LIB.gaia_run.restype = ctypes.c_long
_LIB.gaia_run.argtypes = [
    ctypes.POINTER(ctypes.c_char_p), ctypes.c_int, ctypes.c_long, ctypes.c_int,
    ctypes.c_char_p,
]


def run():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    files = sorted(glob.glob(os.path.join(IN_DIR, "*.gz")))
    n = len(files)
    path_arr = (ctypes.c_char_p * n)(*[p.encode() for p in files])
    written = _LIB.gaia_run(path_arr, n, SLOT, NTHREADS, OUT.encode())
    if written < 0:
        raise RuntimeError("gaia_run failed")
    return written
