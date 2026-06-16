#!/usr/bin/env python3
"""
Preview GWAS SNP file formats to determine correct column names
"""

import pandas as pd
from pathlib import Path
import sys

def preview_files(input_dir):
    
    input_path = Path(input_dir)
    if not input_path.exists():
        print(f"Error: Directory not found: {input_dir}")
        sys.exit(1)
    
    txt_files = list(input_path.glob('*.txt'))
    
    if not txt_files:
        print(f"No .txt files found in {input_dir}")
        sys.exit(1)
    
    print(f"Found {len(txt_files)} .txt files\n")
    print("=" * 80)
    
    for i, file_path in enumerate(txt_files[:5], 1):  # Show first 5 files
        print(f"\nFile {i}: {file_path.name}")
        print("-" * 80)
        
        # Try different delimiters
        for sep_name, sep in [('tab', '\t'), ('space', ' '), ('comma', ',')]:
            try:
                df = pd.read_csv(file_path, sep=sep, nrows=3)
                if len(df.columns) > 1:
                    print(f"\nDelimiter: {sep_name}")
                    print(f"Columns: {', '.join(df.columns.tolist())}")
                    print(f"Shape: {len(df)} rows shown (use head -5 to see more)")
                    print("\nFirst 3 rows:")
                    print(df.to_string(index=False))
                    
                    # Try to identify likely column names
                    cols = df.columns.tolist()
                    print("\nLikely column mappings:")
                    
                    snp_candidates = [c for c in cols if 'snp' in c.lower() or 'rs' in c.lower() or c.upper() in ['SNP', 'RSID', 'RS_ID', 'ID', 'VARIANT_ID']]
                    chr_candidates = [c for c in cols if 'chr' in c.lower() or c.upper() in ['CHR', 'CHROM', 'CHROMOSOME']]
                    pos_candidates = [c for c in cols if 'pos' in c.lower() or 'bp' in c.lower() or c.upper() in ['BP', 'POS', 'POSITION']]
                    p_candidates = [c for c in cols if c.upper() in ['P', 'PVAL', 'P_VALUE', 'PVALUE', 'P.VALUE']]
                    
                    if snp_candidates:
                        print(f"  SNP ID: --snp-col {snp_candidates[0]}")
                    if chr_candidates:
                        print(f"  Chromosome: --chr-col {chr_candidates[0]}")
                    if pos_candidates:
                        print(f"  Position: --pos-col {pos_candidates[0]}")
                    if p_candidates:
                        print(f"  P-value: --p-col {p_candidates[0]}")
                    
                    break
            except Exception as e:
                continue
        
        print("=" * 80)
    
    if len(txt_files) > 5:
        print(f"\n... and {len(txt_files) - 5} more files")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python preview_snp_files.py <directory>")
        sys.exit(1)
    
    preview_files(sys.argv[1])
