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
/// decompression-bound, so the 20 files are gunzipped on a thread pool
/// (Python's gzip releases the GIL, giving true parallelism without the
/// deadlock that an in-process multiprocessing pool hits), largest file first.
/// Parsing the flux arrays and computing min/max/%change is then done
/// vectorized in Polars (native Rust), keeping the scan cost near zero so the
/// wall clock sits close to the raw-decompression floor.
Class Gaia.Benchmark
{

/// Parse the input files and write the challenge CSV.
ClassMethod Process() [ Language = python ]
{
    import sys
    sys.path.insert(0, "/home/irisowner/dev/src")
    import gaia_flux_polars
    n = gaia_flux_polars.run()
    print("{} sources changed by more than 100%".format(n))
}

}
