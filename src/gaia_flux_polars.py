"""Polars-based Gaia DR3 BP/RP flux benchmark (in-process, thread-parallel)."""
import glob
import gzip
import os
from concurrent.futures import ThreadPoolExecutor

import polars as pl

IN_DIR = "/home/irisowner/dev/data/in"
OUT = "/home/irisowner/dev/data/out/challenge_output.csv"


def _decompress(path):
    with open(path, "rb") as f:
        return gzip.decompress(f.read())


def _parse(blob):
    return pl.read_csv(
        blob, comment_prefix="#",
        columns=["source_id", "bp_flux", "rp_flux"],
        schema_overrides={"source_id": pl.Utf8, "bp_flux": pl.Utf8, "rp_flux": pl.Utf8},
    )


def _band(col):
    lst = (pl.col(col).str.strip_chars("[]").str.split(",")
           .list.eval(pl.element().cast(pl.Float64, strict=False)))
    good = lst.list.eval(pl.element().filter(pl.element().is_finite() & (pl.element() > 0)))
    return good.list.min().alias(col + "_min"), good.list.max().alias(col + "_max")


def run():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    # largest file first so the wall-clock-bounding stream starts immediately
    files = sorted(glob.glob(os.path.join(IN_DIR, "*.gz")), key=lambda f: -os.path.getsize(f))
    with ThreadPoolExecutor(os.cpu_count()) as ex:      # gzip releases the GIL -> real parallelism
        blobs = list(ex.map(_decompress, files))
        frames = list(ex.map(_parse, blobs))
    df = pl.concat(frames)
    out = (df.with_columns(*_band("bp_flux"), *_band("rp_flux"))
             .with_columns(
                 bp_pc=(pl.col("bp_flux_max") - pl.col("bp_flux_min")) / pl.col("bp_flux_min") * 100,
                 rp_pc=(pl.col("rp_flux_max") - pl.col("rp_flux_min")) / pl.col("rp_flux_min") * 100)
             .with_columns(percentage_change=pl.max_horizontal("bp_pc", "rp_pc"))
             .filter(pl.col("percentage_change") > 100)
             .sort("percentage_change", descending=True)
             .select("source_id", "bp_flux_min", "bp_flux_max", "rp_flux_min", "rp_flux_max", "percentage_change")
             .rename({"bp_flux_min": "bp_min_flux", "bp_flux_max": "bp_max_flux",
                      "rp_flux_min": "rp_min_flux", "rp_flux_max": "rp_max_flux"}))
    out.write_csv(OUT)
    return out.height
