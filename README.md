# InterSystems Programming Challenge 1: Benchmark

Speed-optimized solution to the challenge: identify Gaia DR3 sources whose **BP or RP flux changed by more than 100%** over the observation period and write the result CSV. Optimized for **shortest execution time** (see the sibling repo `gaia-codegolf` for the minimal-code variant).

## Build & run

```bash
docker-compose up --build -d
docker-compose exec iris iris session iris
USER>do ^RunScript
```

This writes `data/out/challenge_output.csv`.

## How it works

`do ^RunScript` (`src/RunScript.mac`) calls `Gaia.Benchmark.Process()` (`src/Gaia/Benchmark.cls`), which runs the whole job **in-process** in Embedded Python by invoking a compiled C kernel (`src/gaiascan.c` → `gaiascan.so`). Python only passes the input directory and output path; every heavy step happens in C with the GIL released:

- **One fused pass per file, in parallel.** Each OpenMP thread takes one `.csv.gz` file and, in a single pass, `mmap`s it, decompresses it with **libdeflate** (the fastest gzip decoder — the output buffer is sized exactly from the gzip trailer's ISIZE field, so there is one allocation and no resizing), then immediately scans its rows into that file's output slot. There is no separate row-indexing pass and no second parallel region: profiling showed the scan is ~10× cheaper than the decompress, so splitting it across threads bought nothing while the index pass cost more than the scan it was meant to accelerate. For every row the kernel jumps straight to the three needed columns (`source_id`, `bp_flux`, `rp_flux`) with `memchr`, then walks each flux array once with a fast custom float parser, tracking the min and max of the finite, positive values — no intermediate storage of individual fluxes.
- **Largest file first.** The job is bounded by the largest (unsplittable) gzip stream, so the kernel sorts files by size and schedules the biggest first. That starts the longest decompress at t=0 instead of stranding one thread at the end after everything else has finished. With that in place, ~0.75× cores is the fastest thread count — extra threads only add memory-bandwidth contention on the parallel decompress.
- **Compute, filter, sort, write.** `percentage_change = ((max - min) / min) * 100` per band, the larger of BP and RP is kept, rows above 100% are sorted descending, and the CSV is formatted in parallel (each row into its own slot) and written in a single call.

The kernel is compiled at image build with `-O3 -march=native -funroll-loops`; if `gcc`/libdeflate are unavailable it lazily recompiles on first run, and if that fails it falls back to an equivalent pure-Polars engine (`src/gaia_flux_polars.py`) so the run always succeeds.

Output columns: `source_id`, `bp_min_flux`, `bp_max_flux`, `rp_min_flux`, `rp_max_flux`, `percentage_change`.

## Valid Flux Assumptions

The spec says to ignore *"missing, null, NaN, or otherwise invalid"* flux values. This solution treats a flux as valid only if it is **present, numeric, finite, and strictly positive** — i.e. it drops:

- missing / `null` values
- `NaN` / non-finite values
- zero and negative values

If a band has no valid fluxes, its `min`/`max` cells are left empty and only the other band contributes to `percentage_change`.

## Result

**~0.21 seconds** to process all 20 files (1.54 GB uncompressed).

![Benchmark results](benchmark_results.png)
