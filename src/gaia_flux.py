'''Parallel parser + writer for the Gaia DR3 BP/RP flux benchmark.

Run as a standalone script (not imported into the IRIS server process, where a
multiprocessing Pool would deadlock):

    irispython gaia_flux.py <in_dir> <out_csv>

It parses every EpochPhotometry_*.csv.gz in <in_dir> concurrently -- one worker
process per file -- filters to sources whose BP/RP flux %change exceeds 100%,
sorts them, and writes the challenge CSV to <out_csv> directly. Writing the CSV
here (rather than handing rows back to IRIS as JSON) avoids a serialize/readback
round-trip.
'''
import csv
import glob
import gzip
import json
import math
import os
import sys

csv.field_size_limit(10 ** 7)

THRESHOLD = 100  # report sources whose flux changed by more than this percent


def _band(raw):
    '''(min, max, pct_change) over valid (finite, positive) fluxes, or None.'''
    if not raw:
        return None
    try:
        values = json.loads(raw)
    except Exception:
        return None
    kept = [v for v in values
            if isinstance(v, (int, float)) and math.isfinite(v) and v > 0]
    if len(kept) < 2:
        return None
    mn, mx = min(kept), max(kept)
    return mn, mx, ((mx - mn) / mn) * 100.0


def process_file(path):
    '''Rows above the threshold in one file: [(source_id, bp_min, bp_max, rp_min, rp_max, pct), ...].

    A missing band's flux columns are emitted as "" (empty CSV cell).
    '''
    rows = []
    with gzip.open(path, 'rt') as handle:
        reader = csv.DictReader(line for line in handle if not line.startswith('#'))
        for row in reader:
            bp = _band(row.get('bp_flux'))
            rp = _band(row.get('rp_flux'))
            if bp is None and rp is None:
                continue
            pct = max(b[2] for b in (bp, rp) if b)
            if pct <= THRESHOLD:
                continue
            rows.append((
                int(row['source_id']),
                bp[0] if bp else "", bp[1] if bp else "",
                rp[0] if rp else "", rp[1] if rp else "",
                pct,
            ))
    return rows


if __name__ == '__main__':
    from multiprocessing import Pool

    in_dir, out_csv = sys.argv[1], sys.argv[2]
    files = sorted(glob.glob(os.path.join(in_dir, 'EpochPhotometry_*.csv.gz')))
    # One worker per file: with 20 uneven files on ~16 cores, matching the pool
    # size to the file count lets every file start at once (no second wave where
    # a few workers each handle a 2nd file while the rest sit idle).
    with Pool(len(files)) as pool:
        per_file = pool.map(process_file, files)
    rows = sorted((row for chunk in per_file for row in chunk), key=lambda r: -r[-1])
    with open(out_csv, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(["source_id", "bp_min_flux", "bp_max_flux",
                         "rp_min_flux", "rp_max_flux", "percentage_change"])
        writer.writerows(rows)
    print('{} files, {} sources > {}%'.format(len(files), len(rows), THRESHOLD))
