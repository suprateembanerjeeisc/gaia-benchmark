# InterSystems Programming Challenge 1: Benchmark

Speed-optimized solution to the challenge: identify Gaia DR3 sources whose **BP or RP flux changed by more than 100%** over the observation period and write the result CSV. Optimized for **shortest execution time** — sub-second warm (~0.8 s) — by treating the job as what it is: decompression-bound (see the sibling repo `gaia-codegolf` for the minimal-code variant).

## Build & run

```bash
docker-compose up --build -d
docker-compose exec iris iris session iris
USER>do ^RunScript
```

This writes `data/out/challenge_output.csv` (~0.8 s warm).

## How it works

`do ^RunScript` (`src/RunScript.mac`) calls `Gaia.Benchmark.Process()` (`src/Gaia/Benchmark.cls`), which runs `src/gaia_flux_polars.py` **entirely in-process** in Embedded Python — no subprocess, no marshalling of rows back into IRIS. The pipeline:

- **Decompress on a thread pool, largest file first.** Profiling showed the job is decompression-bound, and Python's `gzip` releases the GIL — so the 20 `.csv.gz` files gunzip in genuine parallel across the container's cores. (Threads, not a `multiprocessing` pool: a pool deadlocks when launched inside the IRIS server process; threads do not.) Starting the biggest file first means it — which alone bounds the wall clock — never waits.
- **Parse + compute vectorized in Polars (native Rust).** The raw CSV bytes go straight to Polars, which reads only the three needed columns, splits each `bp_flux` / `rp_flux` array, keeps the valid (finite, positive) values, and computes `min` / `max` and `((max - min) / min) * 100` per band — all as native column operations, keeping the scan cost near zero.
- **Filter, sort, write.** Keep `percentage_change > 100` (the larger of the BP and RP change), sort descending, and write the CSV directly from Polars.

Output columns: `source_id`, `bp_min_flux`, `bp_max_flux`, `rp_min_flux`, `rp_max_flux`, `percentage_change`.

## Valid Flux Assumptions

The spec says to ignore *"missing, null, NaN, or otherwise invalid"* flux values. This solution treats a flux as valid only if it is **present, numeric, finite, and strictly positive** — i.e. it drops:

- missing / `null` values
- `NaN` / non-finite values (Polars `is_finite`)
- zero and negative values

If a band has no valid fluxes, its `min`/`max` cells are left empty and only the
other band contributes to `percentage_change`.

## Verifiable Result

- **~0.8 seconds** warm to process all 20 files (down from ~12 s single-threaded).