#!/usr/bin/env python3
"""
Build wide PRS matrix across all traits.
Matches the wide matrix building logic from the Jupyter notebook.
"""
import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np

def parse_args():
    ap = argparse.ArgumentParser(description="Build wide PRS matrix across traits")
    ap.add_argument("--out", required=True, help="Output matrix file")
    ap.add_argument("--in", dest="inputs", action="append", default=[], 
                    help="Repeatable: --in trait.genomewide.sscore.tsv")
    return ap.parse_args()

def main():
    args = parse_args()
    out = Path(args.out)
    inputs = [
        Path(p) for p in args.inputs 
        if p and Path(p).exists() and Path(p).stat().st_size > 0
    ]
    
    print(f"[build_prs_matrix] Building wide matrix from {len(inputs)} trait files")
    
    if not inputs:
        print(f"[build_prs_matrix] WARNING: No valid input files")
        out.write_text("")
        return 0
    
    wide = None
    used = 0
    skipped = 0
    trait_names = []
    
    for idx, path in enumerate(inputs, 1):
        # Extract trait name - remove .genomewide.sscore.tsv suffix
        trait = path.name.replace(".genomewide.sscore.tsv", "")
        
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception as e:
            print(f"[build_prs_matrix] ERROR reading {path.name}: {e}", file=sys.stderr)
            skipped += 1
            continue
        
        if df.empty:
            print(f"[build_prs_matrix] WARNING: {path.name} is empty")
            skipped += 1
            continue
        
        # Expect column named <trait>_SUM - exactly as in notebook
        col = f"{trait}_SUM"
        if col not in df.columns:
            print(f"[build_prs_matrix] WARNING: {path.name} missing column '{col}'")
            print(f"  Available columns: {list(df.columns)}")
            skipped += 1
            continue
        
        # Verify FID/IID present
        if "FID" not in df.columns or "IID" not in df.columns:
            print(f"[build_prs_matrix] WARNING: {path.name} missing FID/IID")
            skipped += 1
            continue
        
        # Rename <trait>_SUM to just <trait> for compact matrix
        df = df.rename(columns={col: trait})[["FID", "IID", trait]]
        
        # Convert to numeric
        df[trait] = pd.to_numeric(df[trait], errors="coerce")
        
        # Count non-missing values
        n_valid = df[trait].notna().sum()
        print(f"[build_prs_matrix] {idx}/{len(inputs)}: {trait} "
              f"({len(df)} samples, {n_valid} with scores)")
        
        if wide is None:
            wide = df
        else:
            # Outer merge to keep all samples - exactly as in notebook
            wide = wide.merge(df, on=["FID", "IID"], how="outer")
        
        used += 1
        trait_names.append(trait)
    
    print(f"[build_prs_matrix] Summary: used={used}, skipped={skipped}")
    
    if wide is None or wide.empty:
        print(f"[build_prs_matrix] WARNING: No data to write")
        out.write_text("")
        return 0
    
    # Report statistics
    print(f"[build_prs_matrix] Final matrix shape: {wide.shape}")
    print(f"[build_prs_matrix] Traits: {', '.join(trait_names)}")
    
    # Count missing values per trait
    print(f"[build_prs_matrix] Missing value summary:")
    for trait in trait_names:
        if trait in wide.columns:
            n_missing = wide[trait].isna().sum()
            pct_missing = 100 * n_missing / len(wide)
            print(f"  {trait}: {n_missing}/{len(wide)} missing ({pct_missing:.1f}%)")
    
    # Basic statistics for each trait
    print(f"[build_prs_matrix] Per-trait statistics:")
    for trait in trait_names:
        if trait in wide.columns:
            scores = wide[trait].replace([np.inf, -np.inf], np.nan).dropna()
            if len(scores) > 0:
                print(f"  {trait}:")
                print(f"    N: {len(scores)}")
                print(f"    Mean: {scores.mean():.6f}")
                print(f"    Std: {scores.std():.6f}")
                print(f"    Range: [{scores.min():.6f}, {scores.max():.6f}]")
    
    # Write output - exactly as in notebook
    wide.to_csv(out, sep="\t", index=False)
    print(f"[build_prs_matrix] Wrote matrix to {out}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
