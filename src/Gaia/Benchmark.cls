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
/// The 20 files are parsed in parallel (one process per file) by gaia_flux.py,
/// run as a standalone irispython process because a multiprocessing Pool cannot
/// run inside the IRIS server process (it deadlocks). That script also writes
/// the CSV directly, so no rows are marshalled back into IRIS.
Class Gaia.Benchmark
{

/// Spawn the parallel parser, which parses + filters + writes the challenge CSV.
ClassMethod Process() [ Language = python ]
{
    import os, subprocess

    os.makedirs("/home/irisowner/dev/data/out", exist_ok=True)
    subprocess.run(
        ["irispython",
         "/home/irisowner/dev/src/gaia_flux.py",
         "/home/irisowner/dev/data/in",
         "/home/irisowner/dev/data/out/challenge_output.csv"],
        check=True,
    )
}

}
