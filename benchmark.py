#!/usr/bin/env python3
"""Repeatedly run `do ^RunScript` in the benchmark container and record timings.

Runs the challenge N times (default 1000), writing every run's elapsed time to
benchmark.csv and printing a summary distribution at the end.

Usage:
    python benchmark.py            # 1000 runs
    python benchmark.py 200        # 200 runs
    python benchmark.py --container challenge-submission-iris-1
"""
import argparse
import csv
import re
import statistics
import subprocess
import sys
import time

DEFAULT_RUNS = 1000
DEFAULT_CONTAINER = "challenge-submission-iris-1"
OUTPUT_CSV = "benchmark.csv"
ELAPSED_RE = re.compile(r"Elapsed time:\s*([\d.]+)")


def one_run(container):
    """Run `do ^RunScript` once; return the routine's reported elapsed seconds."""
    proc = subprocess.run(
        ["docker", "exec", "-i", container, "iris", "session", "iris"],
        input="do ^RunScript\nhalt\n",
        capture_output=True, text=True, timeout=300,
    )
    match = ELAPSED_RE.search(proc.stdout)
    if not match:
        raise RuntimeError(f"no elapsed time in output:\n{proc.stdout[-500:]}")
    return float(match.group(1))


def summarize(times):
    ordered = sorted(times)

    def pct(p):
        k = (len(ordered) - 1) * p
        lo = int(k)
        hi = min(lo + 1, len(ordered) - 1)
        return ordered[lo] + (ordered[hi] - ordered[lo]) * (k - lo)

    print(f"\nruns    = {len(times)}")
    print(f"min     = {min(times):.3f}s")
    print(f"max     = {max(times):.3f}s")
    print(f"mean    = {statistics.mean(times):.3f}s")
    print(f"median  = {statistics.median(times):.3f}s")
    print(f"stdev   = {statistics.pstdev(times):.3f}s")
    print(f"p10     = {pct(0.10):.3f}s")
    print(f"p90     = {pct(0.90):.3f}s")
    print(f"p95     = {pct(0.95):.3f}s")
    print(f"p99     = {pct(0.99):.3f}s")

    print("\nhistogram (50 ms bins):")
    bins = {}
    for x in times:
        b = round(x // 0.05 * 0.05, 2)
        bins[b] = bins.get(b, 0) + 1
    peak = max(bins.values())
    for b in sorted(bins):
        bar = "#" * max(1, round(bins[b] / peak * 40))
        print(f"  {b:.2f}-{b + 0.05:.2f}s | {bar} ({bins[b]})")


def main():
    ap = argparse.ArgumentParser(description="Benchmark do ^RunScript over many runs.")
    ap.add_argument("runs", nargs="?", type=int, default=DEFAULT_RUNS,
                    help=f"number of runs (default {DEFAULT_RUNS})")
    ap.add_argument("--container", default=DEFAULT_CONTAINER,
                    help=f"IRIS container name (default {DEFAULT_CONTAINER})")
    ap.add_argument("--csv", default=OUTPUT_CSV, help=f"output CSV (default {OUTPUT_CSV})")
    args = ap.parse_args()

    times = []
    started = time.time()
    with open(args.csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["run", "elapsed_seconds"])
        for i in range(1, args.runs + 1):
            try:
                elapsed = one_run(args.container)
            except Exception as exc:
                print(f"\nrun {i} failed: {exc}", file=sys.stderr)
                continue
            times.append(elapsed)
            writer.writerow([i, elapsed])
            fh.flush()  # persist as we go, so a long run is never lost
            # single-line progress: running distribution, not the noisy last value
            print(f"\r{i}/{args.runs}  min={min(times):.3f}  max={max(times):.3f}  "
                  f"mean={statistics.mean(times):.3f}  median={statistics.median(times):.3f}",
                  end="", flush=True)

    print(f"\n\nwrote {len(times)} rows to {args.csv} "
          f"({time.time() - started:.0f}s wall)")
    if times:
        summarize(times)


if __name__ == "__main__":
    main()
