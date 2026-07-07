"""Polars-based Gaia DR3 BP/RP flux benchmark (in-process, thread-parallel).

Decompression uses libdeflate (the fastest gzip decoder) via ctypes, which
releases the GIL on the foreign call so the 20 files gunzip in true parallel
across threads. If libdeflate can't be loaded, it falls back to the stdlib
``gzip`` module (also GIL-releasing, just slower). Parsing and the min/max /
percentage-change computation are done vectorized in Polars.
"""
import ctypes
import ctypes.util
import glob
import gzip
import os
import struct
from concurrent.futures import ThreadPoolExecutor

import polars as pl

IN_DIR = "/home/irisowner/dev/data/in"
OUT = "/home/irisowner/dev/data/out/challenge_output.csv"


def _load_libdeflate():
    """Return a configured libdeflate handle, or None if unavailable."""
    for name in ("deflate", "libdeflate.so.0", "libdeflate.so", ctypes.util.find_library("deflate")):
        if not name:
            continue
        try:
            lib = ctypes.CDLL(name)
            lib.libdeflate_alloc_decompressor.restype = ctypes.c_void_p
            lib.libdeflate_alloc_decompressor.argtypes = []
            lib.libdeflate_free_decompressor.argtypes = [ctypes.c_void_p]
            lib.libdeflate_gzip_decompress.argtypes = [
                ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t,
                ctypes.c_char_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
            lib.libdeflate_gzip_decompress.restype = ctypes.c_int
            return lib
        except (OSError, AttributeError):
            continue
    return None


_LIBDEFLATE = _load_libdeflate()


def _decompress(path):
    with open(path, "rb") as f:
        blob = f.read()
    if _LIBDEFLATE is None:
        return gzip.decompress(blob)
    # gzip trailer's ISIZE (last 4 bytes) gives the exact output size -> one alloc.
    size = struct.unpack("<I", blob[-4:])[0]
    dec = _LIBDEFLATE.libdeflate_alloc_decompressor()
    try:
        out = ctypes.create_string_buffer(size)
        actual = ctypes.c_size_t(0)
        rc = _LIBDEFLATE.libdeflate_gzip_decompress(
            dec, blob, len(blob), out, size, ctypes.byref(actual))
        if rc != 0:
            return gzip.decompress(blob)  # fall back on any decode error
        return out.raw[:actual.value]
    finally:
        _LIBDEFLATE.libdeflate_free_decompressor(dec)


def _band(col):
    lst = (pl.col(col).str.strip_chars("[]").str.split(",")
           .list.eval(pl.element().cast(pl.Float64, strict=False)))
    good = lst.list.eval(pl.element().filter(pl.element().is_finite() & (pl.element() > 0)))
    return good.list.min().alias(col + "_min"), good.list.max().alias(col + "_max")


def _process(path):
    """Decompress + parse + reduce one file to its qualifying rows.

    Fused per-file so each file's decompress -> parse -> min/max -> filter runs
    as one unit in the thread pool. Parsing of already-decompressed files thus
    overlaps decompression of the still-running (largest) ones, instead of
    waiting on a barrier for every file to finish decompressing first.
    """
    df = pl.read_csv(
        _decompress(path), comment_prefix="#",
        columns=["source_id", "bp_flux", "rp_flux"],
        schema_overrides={"source_id": pl.Utf8, "bp_flux": pl.Utf8, "rp_flux": pl.Utf8},
        rechunk=False,
    )
    return (df.with_columns(*_band("bp_flux"), *_band("rp_flux"))
              .with_columns(
                  bp_pc=(pl.col("bp_flux_max") - pl.col("bp_flux_min")) / pl.col("bp_flux_min") * 100,
                  rp_pc=(pl.col("rp_flux_max") - pl.col("rp_flux_min")) / pl.col("rp_flux_min") * 100)
              .with_columns(percentage_change=pl.max_horizontal("bp_pc", "rp_pc"))
              .filter(pl.col("percentage_change") > 100)
              .select("source_id",
                      pl.col("bp_flux_min").alias("bp_min_flux"),
                      pl.col("bp_flux_max").alias("bp_max_flux"),
                      pl.col("rp_flux_min").alias("rp_min_flux"),
                      pl.col("rp_flux_max").alias("rp_max_flux"),
                      "percentage_change"))


def run():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    # largest file first so the wall-clock-bounding stream starts immediately
    files = sorted(glob.glob(os.path.join(IN_DIR, "*.gz")), key=lambda f: -os.path.getsize(f))
    with ThreadPoolExecutor(os.cpu_count()) as ex:  # libdeflate/gzip release the GIL -> real parallelism
        parts = list(ex.map(_process, files))
    out = pl.concat(parts).sort("percentage_change", descending=True)
    out.write_csv(OUT)
    return out.height
