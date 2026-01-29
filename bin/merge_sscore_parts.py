#!/usr/bin/env python3
"""
Merge per-chromosome .sscore files into genome-wide PRS.
Matches the aggregation logic from the Jupyter notebook.
"""
import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np

def parse_args():
    ap = argparse.ArgumentParser(description="Merge per-chromosome sscore files")
    ap.add_argument("--trait", required=True, help="Trait name")
    ap.add_argument("--out", required=True, help="Output file path")
    ap.add_argument("--part", action="append", default=[], 
                    help="Repeatable: --part file.sscore")
    return ap.parse_args()

def main():
    args = parse_args()
    trait = args.trait
    out = Path(args.out)
    parts = [Path(p) for p in args.part if p and Path(p).exists()]
    
    print(f"[merge_sscore_parts] Trait: {trait}")
    print(f"[merge_sscore_parts] Found {len(parts)} part files")
    
    if not parts:
        print(f"[merge_sscore_parts] WARNING: No valid part files")
        out.write_text("")
        return 0
    
    merged = None
    used = 0
    skipped_empty = 0
    skipped_nocols = 0
    skipped_noscore = 0
    
    for idx, part_path in enumerate(parts, 1):
        # Check file size
        if part_path.stat().st_size == 0:
            skipped_empty += 1
            continue
        
        try:
            df = pd.read_csv(part_path, sep="\t")
        except Exception as e:
            print(f"[merge_sscore_parts] ERROR reading {part_path.name}: {e}", 
                  file=sys.stderr)
            continue
        
        if df.empty:
            skipped_empty += 1
            continue
        
        # Normalize FID/IID columns - exactly as in notebook
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
            print(f"[merge_sscore_parts] WARNING: {part_path.name} missing FID/IID")
            skipped_nocols += 1
            continue
        
        # Find PRS column - prefer SCORE_SUM, accept SCORE1_SUM
        prs_col = None
        for c in df.columns:
            if c.strip().upper() == "SCORE_SUM":
                prs_col = c
                break
        if prs_col is None:
            for c in df.columns:
                if c.strip().upper() == "SCORE1_SUM":
                    prs_col = c
                    break
        
        if prs_col is None:
            print(f"[merge_sscore_parts] WARNING: {part_path.name} has no SCORE_SUM column")
            skipped_noscore += 1
            continue
        
        # Extract chromosome label from filename - exactly as in notebook
        fn = part_path.name
        chr_label = "NA"
        if ".chr" in fn and ".scored" in fn:
            chr_label = fn.split(".chr")[1].split(".scored")[0]
        
        # Keep only FID, IID, and score
        df_small = df[["FID", "IID", prs_col]].copy()
        df_small = df_small.rename(columns={prs_col: f"SCORE_SUM_chr{chr_label}"})
        
        # Convert to numeric and handle missing values
        df_small[f"SCORE_SUM_chr{chr_label}"] = pd.to_numeric(
            df_small[f"SCORE_SUM_chr{chr_label}"], errors="coerce"
        )
        
        if merged is None:
            merged = df_small
        else:
            # Outer merge to keep all samples - exactly as in notebook
            merged = merged.merge(df_small, on=["FID", "IID"], how="outer")
        
        used += 1
        
        if idx % 5 == 0 or idx == len(parts):
            print(f"[merge_sscore_parts] Merged {idx}/{len(parts)} parts, "
                  f"current shape: {merged.shape if merged is not None else (0,0)}")
    
    print(f"[merge_sscore_parts] Summary: used={used}, skipped_empty={skipped_empty}, "
          f"skipped_nocols={skipped_nocols}, skipped_noscore={skipped_noscore}")
    
    if merged is None or merged.empty:
        print(f"[merge_sscore_parts] WARNING: No usable parts, writing empty file")
        out.write_text("")
        return 0
    
    # Sum across chromosomes - exactly as in notebook
    score_cols = [c for c in merged.columns if c.startswith("SCORE_SUM_chr")]
    if not score_cols:
        print(f"[merge_sscore_parts] ERROR: No SCORE_SUM_chr* columns after merge")
        out.write_text("")
        return 1
    
    print(f"[merge_sscore_parts] Summing {len(score_cols)} chromosome scores")
    
    # Sum with skipna=True - exactly as in notebook
    merged[f"{trait}_SUM"] = merged[score_cols].sum(axis=1, skipna=True)
    
    # Check for samples with all NaN (no scores)
    all_nan = merged[score_cols].isna().all(axis=1).sum()
    if all_nan > 0:
        print(f"[merge_sscore_parts] WARNING: {all_nan} samples have no scores across all chromosomes")
    
    # Statistics on final scores
    final_scores = merged[f"{trait}_SUM"].replace([np.inf, -np.inf], np.nan).dropna()
    if len(final_scores) > 0:
        print(f"[merge_sscore_parts] Final PRS statistics:")
        print(f"  N samples: {len(final_scores)}")
        print(f"  Mean: {final_scores.mean():.6f}")
        print(f"  Std: {final_scores.std():.6f}")
        print(f"  Min: {final_scores.min():.6f}")
        print(f"  Max: {final_scores.max():.6f}")
    
    # Write output - exactly as in notebook
    merged[["FID", "IID", f"{trait}_SUM"]].to_csv(out, sep="\t", index=False)
    print(f"[merge_sscore_parts] Wrote {len(merged)} samples to {out}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

