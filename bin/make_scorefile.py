#!/usr/bin/env python3
"""
Create PLINK2 score file from GWAS weights + PVAR.
Implements your notebook logic:
  - merge by (chr,pos)
  - direct vs swapped allele match
  - optional strand-flip
Outputs: score.tsv with columns: ID ALLELE SCORE
"""

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--weights", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--pvar", required=True)
    ap.add_argument("--out_score", required=True)
    ap.add_argument(
        "--score_name",
        default="SCORE",
        help="Name for the coefficient column in the output score file header."
    )
    ap.add_argument("--allow_strand_flip", action="store_true")
    return ap.parse_args()

def norm_chr(x: str) -> str:
    return str(x).strip().replace("chr", "")

def main():
    args = parse_args()
    chrom = str(args.chrom)

    wpath = Path(args.weights)
    pvar_path = Path(args.pvar)
    out_score = Path(args.out_score)

    print(f"[make_scorefile] weights={wpath}")
    print(f"[make_scorefile] chrom={chrom}")
    print(f"[make_scorefile] pvar={pvar_path}")
    print(f"[make_scorefile] out={out_score}")
    print(f"[make_scorefile] allow_strand_flip={args.allow_strand_flip}")
    print(f"[make_scorefile] score_name={args.score_name}")

    # Read weights
    w = pd.read_csv(wpath, sep="\t", compression="infer")
    req = {"chromosome", "base_pair_location", "effect_allele", "other_allele", "beta"}
    missing = sorted(list(req - set(w.columns)))
    if missing:
        print(f"ERROR: weights missing columns: {missing}", file=sys.stderr)
        out_score.write_text("")
        return 1

    w["chromosome"] = w["chromosome"].astype(str).map(norm_chr)
    w["base_pair_location"] = w["base_pair_location"].astype(int)
    w["other_allele"] = w["other_allele"].astype(str).str.upper().str.strip()
    w["effect_allele"] = w["effect_allele"].astype(str).str.upper().str.strip()
    w["beta"] = pd.to_numeric(w["beta"], errors="coerce")

    w = w[w["chromosome"] == chrom].copy()
    w = w.dropna(subset=["beta"])
    if w.empty:
        print("[make_scorefile] no weights for this chrom -> empty output")
        out_score.write_text("")
        return 0

    # Read PVAR (skip header lines beginning with '#')
    pvar = pd.read_csv(
        pvar_path,
        sep="\t",
        comment="#",
        header=None,
        names=["chr", "pos", "id", "ref", "alt"],
        usecols=[0, 1, 2, 3, 4],
    )
    pvar["chr"] = pvar["chr"].astype(str).map(norm_chr)
    pvar["pos"] = pvar["pos"].astype(int)
    pvar["ref"] = pvar["ref"].astype(str).str.upper()
    pvar["alt"] = pvar["alt"].astype(str).str.upper()

    # Keep biallelic only (avoid ALT like A,C)
    pvar = pvar[~pvar["alt"].str.contains(",", regex=False)].copy()

    # Prefilter by positions present in weights
    pos_set = set(w["base_pair_location"].tolist())
    pvar = pvar[pvar["pos"].isin(pos_set)].copy()
    if pvar.empty:
        print("[make_scorefile] no PVAR positions match weights -> empty output")
        out_score.write_text("")
        return 0

    merged = w.merge(
        pvar[["chr", "pos", "id", "ref", "alt"]],
        left_on=["chromosome", "base_pair_location"],
        right_on=["chr", "pos"],
        how="inner",
    )
    if merged.empty:
        print("[make_scorefile] no overlaps after merge -> empty output")
        out_score.write_text("")
        return 0

    direct_ct = swapped_ct = flip_ct = dropped_ct = 0
    rows = []

    for r in merged.itertuples(index=False):
        w_other = r.other_allele
        w_eff = r.effect_allele
        beta = float(r.beta)

        pid = r.id
        pref = r.ref
        palt = r.alt

        matched = False

        # DIRECT
        if w_other == pref and w_eff == palt:
            rows.append((pid, palt, beta))
            direct_ct += 1
            matched = True

        # SWAPPED
        elif w_other == palt and w_eff == pref:
            rows.append((pid, pref, beta))
            swapped_ct += 1
            matched = True

        # STRANDFLIP
        elif args.allow_strand_flip:
            wo_c = COMPLEMENT.get(w_other)
            we_c = COMPLEMENT.get(w_eff)
            if wo_c and we_c:
                if wo_c == pref and we_c == palt:
                    rows.append((pid, palt, beta))
                    flip_ct += 1
                    matched = True
                elif wo_c == palt and we_c == pref:
                    rows.append((pid, pref, beta))
                    flip_ct += 1
                    matched = True

        if not matched:
            dropped_ct += 1

    print(f"[make_scorefile] matched: direct={direct_ct}, swapped={swapped_ct}, flip={flip_ct}, dropped={dropped_ct}")

    if not rows:
        out_score.write_text("")
        return 0

    # Use a per-trait score name so PLINK2 can carry it through to .sscore column names.
    score_name = str(args.score_name).strip() or "SCORE"
    score_df = pd.DataFrame(rows, columns=["ID", "ALLELE", score_name])

    # Deduplicate variant IDs by |beta| (same as your notebook tie-breaker)
    if score_df["ID"].duplicated().any():
        score_df["absS"] = score_df[score_name].abs()
        score_df = (
            score_df.sort_values(["ID", "absS"], ascending=[True, False])
                    .drop_duplicates("ID", keep="first")
                    .drop(columns="absS")
        )

    score_df.to_csv(out_score, sep="\t", index=False)
    print(f"[make_scorefile] wrote {len(score_df)} score variants -> {out_score}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())


