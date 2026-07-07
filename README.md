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

- **Phase 1 — parallel decompress.** All 20 `.csv.gz` files are decompressed with **libdeflate** (the fastest gzip decoder), one file per OpenMP thread. The output buffer for each file is sized exactly from the gzip trailer's ISIZE field, so there is a single allocation and no resizing.
- **Phase 2 — parallel single-pass scan.** Each decompressed buffer is split into row-chunks, and *all* chunks across *all* files are scanned in one parallel loop — so the largest file's scan spreads across threads instead of becoming a straggler. For every row the kernel jumps directly to the three needed columns (`source_id`, `bp_flux`, `rp_flux`) with `memchr`, then walks each flux array once with a fast custom float parser, tracking the min and max of the finite, positive values. No intermediate storage of individual fluxes.
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

**~0.29 seconds** to process all 20 files.

![Benchmark results](benchmark_results.png)
