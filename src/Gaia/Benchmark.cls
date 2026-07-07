/// Gaia DR3 BP/RP flux-variability benchmark.
///
/// Process() reads the 20 EpochPhotometry_*.csv.gz files shipped in ./data/in,
/// computes each source's BP/RP flux percentage change, and writes the objects
/// whose change exceeds 100% to ./data/out/challenge_output.csv.
///
/// Percentage change for a band = ((max_flux - min_flux) / min_flux) * 100,
/// over the valid (finite, positive) flux values; percentage_change is the
/// larger of the BP and RP values.
///
/// The work runs entirely in-process in Embedded Python. The job is
/// decompression-bound, so both engines gunzip all 20 files in parallel with
/// the GIL released:
///   - gaia_fast: a C kernel (gaiascan.so) that decompresses with libdeflate
///     and scans each file's BP/RP arrays in a single pass (OpenMP over files),
///     returning only per-source min/max scalars — the fastest path.
///   - gaia_flux_polars: pure Python/Polars (libdeflate via ctypes + vectorized
///     parse) — the portable fallback if the C kernel can't be built/loaded.
/// Both share the same percentage-change / filter / sort / write tail, so the
/// output is identical.
Class Gaia.Benchmark
{

/// Parse the input files and write the challenge CSV.
ClassMethod Process() [ Language = python ]
{
    import sys
    sys.path.insert(0, "/home/irisowner/dev/src")
    try:
        import gaia_fast as engine
    except Exception as exc:
        print("fast engine unavailable ({}); using Polars engine".format(exc))
        import gaia_flux_polars as engine
    n = engine.run()
    print("{} sources changed by more than 100%".format(n))
}

}
