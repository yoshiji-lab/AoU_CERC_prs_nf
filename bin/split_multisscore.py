#!/usr/bin/env python3
"""Split a multi-score PLINK2 .sscore file into per-trait .sscore parts.

We use PLINK2 --score-list to produce a single .sscore with many <score>_SUM
columns. Downstream, merge_sscore_parts.py expects per-trait per-chrom files
with a SCORE_SUM column.

This script:
  - Reads the multi-score .sscore
  - For each <name>_SUM column, writes <name>.chr<chrom>.scored.sscore
    containing FID, IID, SCORE_SUM.
  - Optionally creates empty .sscore files for traits that were not scored
    on this chromosome (e.g., no overlap).
"""

import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sscore", required=True, help="Multi-score .sscore from PLINK2")
    ap.add_argument("--chrom", required=True, help="Chromosome label for filenames")
    ap.add_argument(
        "--expected_traits",
        default=None,
        help="Optional: text file with one trait name per line. Empty .sscore files will be created for traits missing from the multi-score output.",
    )
    return ap.parse_args()


def main():
    args = parse_args()
    sscore_path = Path(args.sscore)
    chrom = str(args.chrom)

    expected = set()
    if args.expected_traits:
        for line in Path(args.expected_traits).read_text().splitlines():
            t = line.strip()
            if t:
                expected.add(t)

    if (not sscore_path.exists()) or sscore_path.stat().st_size == 0:
        # Nothing scored; just emit empty files if expected list provided.
        for t in sorted(expected):
            Path(f"{t}.chr{chrom}.scored.sscore").write_text("")
        return 0

    df = pd.read_csv(sscore_path, sep="\t")
    if df.empty:
        for t in sorted(expected):
            Path(f"{t}.chr{chrom}.scored.sscore").write_text("")
        return 0

    # Normalize FID/IID like merge_sscore_parts.py does.
    rename_map = {}
    for c in df.columns:
        c_norm = c.strip().lower()
        if c_norm == "#fid" or c_norm == "fid":
            rename_map[c] = "FID"
        elif c_norm == "iid":
            rename_map[c] = "IID"
    if rename_map:
        df = df.rename(columns=rename_map)

    if "FID" not in df.columns or "IID" not in df.columns:
        raise SystemExit(f"ERROR: {sscore_path} missing FID/IID columns")

    # Identify score sum columns.
    sum_cols = [c for c in df.columns if c.endswith("_SUM")]
    produced = set()

    for c in sum_cols:
        trait = c[: -len("_SUM")]
        produced.add(trait)
        out = Path(f"{trait}.chr{chrom}.scored.sscore")
        out_df = df[["FID", "IID", c]].copy()
        out_df = out_df.rename(columns={c: "SCORE_SUM"})
        out_df.to_csv(out, sep="\t", index=False)

    # Create empty files for missing expected traits.
    for t in sorted(expected - produced):
        Path(f"{t}.chr{chrom}.scored.sscore").write_text("")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

