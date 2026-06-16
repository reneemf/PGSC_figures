#!/usr/bin/env python3
"""
GWAS SNP Annotation Pipeline
Annotates SNPs with gene locations (±10kb) and ClinVar clinical significance
Build: GRCh37/hg19
"""

import requests
import pandas as pd
import time
import argparse
from pathlib import Path
import json
from typing import Dict, List, Optional
import sys

class SNPAnnotator:
    def __init__(self, flank_distance=10000, build='hg19'):
        self.flank_distance = flank_distance
        self.build = build
        self.session = requests.Session()
        
    def parse_snp_files(self, input_dir: Path, p_threshold: float, 
                       snp_col='SNP', chr_col='CHR', pos_col='BP', p_col='P',
                       file_patterns=None) -> pd.DataFrame:
        """Parse multiple SNP files and apply Bonferroni correction
        
        Args:
            file_patterns: List of file suffixes to include (e.g., ['_age.txt', '_sex.txt'])
                          If None, processes all .txt files
        """
        all_snps = []
        
        # Get all .txt files
        all_files = list(input_dir.glob('*.txt'))
        
        # Filter by patterns if specified
        if file_patterns:
            files_to_process = [
                f for f in all_files 
                if any(f.name.endswith(pattern) for pattern in file_patterns)
            ]
            print(f"Found {len(files_to_process)} files matching patterns: {', '.join(file_patterns)}")
            print(f"(Skipping {len(all_files) - len(files_to_process)} files)\n")
        else:
            files_to_process = all_files
        
        for file_path in files_to_process:
            print(f"Processing {file_path.name}...")
            try:
                # Try common delimiters - use delim_whitespace for space-separated files
                df = None
                for sep in ['\t', ',']:
                    try:
                        df = pd.read_csv(file_path, sep=sep)
                        if len(df.columns) > 1:
                            break
                    except:
                        continue
                
                # If tab/comma didn't work, try whitespace (handles multiple spaces)
                if df is None or len(df.columns) <= 1:
                    try:
                        df = pd.read_csv(file_path, delim_whitespace=True)
                    except:
                        df = pd.read_csv(file_path, sep=r'\s+', engine='python')
                
                # Check if required columns exist
                if p_col not in df.columns:
                    print(f"  Warning: No P-value column '{p_col}' in {file_path.name}")
                    print(f"  Available columns: {', '.join(df.columns.tolist())}")
                    print(f"  Skipping...")
                    continue
                
                if snp_col not in df.columns:
                    print(f"  Warning: No SNP column '{snp_col}' in {file_path.name}")
                    print(f"  Available columns: {', '.join(df.columns.tolist())}")
                    print(f"  Skipping...")
                    continue
                
                # Apply threshold
                df_sig = df[df[p_col] < p_threshold].copy()
                df_sig['source_file'] = file_path.name
                
                # Add placeholder columns if CHR/BP not present (will be looked up via rsID)
                if chr_col not in df_sig.columns:
                    df_sig[chr_col] = None
                if pos_col not in df_sig.columns:
                    df_sig[pos_col] = None
                
                print(f"  Found {len(df_sig)} SNPs passing threshold from {len(df)} total")
                all_snps.append(df_sig)
                
            except Exception as e:
                print(f"  Error processing {file_path.name}: {e}")
                continue
        
        if not all_snps:
            raise ValueError("No SNPs found passing threshold in any file")
        
        combined = pd.concat(all_snps, ignore_index=True)
        print(f"\nTotal significant SNPs: {len(combined)}")
        
        return combined
    
    def get_position_from_rsid(self, rsid: str) -> tuple:
        """Get chromosome and position from rsID using Ensembl API (GRCh37)"""
        if not rsid or not str(rsid).startswith('rs'):
            return None, None
        
        url = f"https://grch37.rest.ensembl.org/variation/human/{rsid}"
        
        try:
            response = self.session.get(url, headers={"Content-Type": "application/json"})
            response.raise_for_status()
            data = response.json()
            
            if 'mappings' in data and len(data['mappings']) > 0:
                # Get the first mapping (primary assembly)
                mapping = data['mappings'][0]
                chrom = mapping.get('seq_region_name', '').replace('chr', '')
                pos = mapping.get('start', None)
                
                # Some variants have multiple positions; prefer non-MT, non-patch chromosomes
                for m in data['mappings']:
                    seq_name = m.get('seq_region_name', '')
                    if seq_name.isdigit() or seq_name in ['X', 'Y']:
                        chrom = seq_name
                        pos = m.get('start', None)
                        break
                
                return chrom, pos
            
            return None, None
            
        except Exception as e:
            print(f"  Error fetching position for {rsid}: {e}")
            return None, None
        
        time.sleep(0.1)  # Rate limiting
    
    def annotate_genes_ensembl(self, chrom: str, pos: int) -> List[Dict]:
        """Get genes within flanking distance using Ensembl API (GRCh37)"""
        chrom = str(chrom).replace('chr', '')
        start = max(1, pos - self.flank_distance)
        end = pos + self.flank_distance
        
        url = f"https://grch37.rest.ensembl.org/overlap/region/human/{chrom}:{start}-{end}?feature=gene"
        
        try:
            response = self.session.get(url, headers={"Content-Type": "application/json"})
            response.raise_for_status()
            genes = response.json()
            
            results = []
            for gene in genes:
                gene_start = gene.get('start', 0)
                gene_end = gene.get('end', 0)
                
                # Calculate distance
                if pos >= gene_start and pos <= gene_end:
                    distance = 0
                    location = "within"
                elif pos < gene_start:
                    distance = gene_start - pos
                    location = "upstream"
                else:
                    distance = pos - gene_end
                    location = "downstream"
                
                results.append({
                    'gene_name': gene.get('external_name', 'Unknown'),
                    'gene_id': gene.get('id', ''),
                    'biotype': gene.get('biotype', ''),
                    'distance': distance,
                    'location': location,
                    'gene_start': gene_start,
                    'gene_end': gene_end
                })
            
            return sorted(results, key=lambda x: abs(x['distance']))
            
        except Exception as e:
            print(f"  Error annotating genes for {chrom}:{pos}: {e}")
            return []
        
        time.sleep(0.1)  # Rate limiting
    
    def query_clinvar(self, rsid: str = None, chrom: str = None, pos: int = None) -> Dict:
        """Query ClinVar using MyVariant.info API"""
        if rsid:
            query = rsid
        elif chrom and pos:
            chrom = str(chrom).replace('chr', '')
            query = f"chr{chrom}:g.{pos}"
        else:
            return {'clinvar_status': 'not_queried'}
        
        url = f"https://myvariant.info/v1/variant/{query}"
        params = {
            'fields': 'clinvar.rcv.clinical_significance,clinvar.rcv.conditions.name,clinvar.variant_id',
            'assembly': 'hg19'
        }
        
        try:
            response = self.session.get(url, params=params)
            response.raise_for_status()
            data = response.json()
            
            if 'clinvar' in data:
                clinvar = data['clinvar']
                
                # Extract clinical significance
                significances = []
                conditions = []
                
                if 'rcv' in clinvar:
                    rcv_list = clinvar['rcv'] if isinstance(clinvar['rcv'], list) else [clinvar['rcv']]
                    
                    for rcv in rcv_list:
                        if 'clinical_significance' in rcv:
                            sig = rcv['clinical_significance']
                            if isinstance(sig, list):
                                significances.extend(sig)
                            else:
                                significances.append(sig)
                        
                        if 'conditions' in rcv:
                            cond = rcv['conditions']
                            if isinstance(cond, list):
                                for c in cond:
                                    if 'name' in c:
                                        conditions.append(c['name'])
                            elif isinstance(cond, dict) and 'name' in cond:
                                conditions.append(cond['name'])
                
                return {
                    'in_clinvar': True,
                    'clinvar_id': clinvar.get('variant_id', ''),
                    'clinical_significance': '|'.join(set(significances)) if significances else 'Not provided',
                    'conditions': '|'.join(set(conditions)) if conditions else 'Not provided'
                }
            
            return {'in_clinvar': False}
            
        except Exception as e:
            print(f"  Error querying ClinVar for {query}: {e}")
            return {'in_clinvar': False, 'error': str(e)}
        
        time.sleep(0.1)  # Rate limiting
    
    def annotate_snp(self, row: pd.Series, snp_col='SNP', chr_col='CHR', pos_col='BP') -> Dict:
        """Annotate a single SNP with genes and ClinVar info"""
        rsid = row.get(snp_col, None)
        chrom = row.get(chr_col, None)
        pos = row.get(pos_col, None)
        
        # If chr/pos missing but we have rsID, look them up
        if (pd.isna(chrom) or pd.isna(pos)) and rsid and str(rsid).startswith('rs'):
            chrom, pos = self.get_position_from_rsid(rsid)
        
        result = {
            'rsid': rsid,
            'chr': chrom,
            'pos': pos,
        }
        
        # Gene annotation (requires chr and pos)
        if chrom and pos and not pd.isna(chrom) and not pd.isna(pos):
            genes = self.annotate_genes_ensembl(chrom, pos)
            if genes:
                # Get closest gene
                closest = genes[0]
                result['nearest_gene'] = closest['gene_name']
                result['gene_distance'] = closest['distance']
                result['gene_location'] = closest['location']
                result['gene_biotype'] = closest['biotype']
                result['all_nearby_genes'] = '|'.join([g['gene_name'] for g in genes])
            else:
                result['nearest_gene'] = 'None within 10kb'
                result['gene_distance'] = 'N/A'
                result['gene_location'] = 'N/A'
                result['gene_biotype'] = 'N/A'
                result['all_nearby_genes'] = ''
        else:
            result['nearest_gene'] = 'Position unknown'
            result['gene_distance'] = 'N/A'
            result['gene_location'] = 'N/A'
            result['gene_biotype'] = 'N/A'
            result['all_nearby_genes'] = ''
        
        # ClinVar annotation
        clinvar_data = self.query_clinvar(rsid=rsid, chrom=chrom, pos=pos)
        result.update(clinvar_data)
        
        return result
    
    def process_snps(self, snp_df: pd.DataFrame, snp_col='SNP', chr_col='CHR',
                     pos_col='BP', output_file='annotated_snps.csv') -> pd.DataFrame:
        """Process all SNPs and create annotated output"""
        results = []
        total = len(snp_df)
        
        print(f"\nAnnotating {total} SNPs...")
        print("This may take a while due to API rate limiting (~3-5 SNPs per second)")
        
        for idx, row in snp_df.iterrows():
            if idx % 100 == 0:
                print(f"  Progress: {idx}/{total} ({100*idx/total:.1f}%)")
            
            annotation = self.annotate_snp(row, snp_col, chr_col, pos_col)
            
            # Combine original data with annotations
            result_row = row.to_dict()
            result_row.update(annotation)
            results.append(result_row)
        
        print(f"  Completed: {total}/{total} (100%)")
        
        # Create output dataframe
        output_df = pd.DataFrame(results)
        
        # Reorder columns for readability
        priority_cols = ['rsid', 'chr', 'pos', 'nearest_gene', 'gene_distance', 
                        'gene_location', 'in_clinvar', 'clinical_significance', 
                        'conditions', 'source_file']
        other_cols = [c for c in output_df.columns if c not in priority_cols]
        output_df = output_df[[c for c in priority_cols if c in output_df.columns] + other_cols]
        
        # Save results
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        output_df.to_csv(output_file, index=False)
        print(f"\nResults saved to: {output_file}")
        
        return output_df


def main():
    parser = argparse.ArgumentParser(
        description='Annotate GWAS SNPs with gene locations and ClinVar status',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Use Bonferroni threshold (5e-8) with default column names
  python gwas_snp_annotator.py -i /path/to/snps/ -o annotated_results.csv

  # Only process files ending with _age.txt, _sex.txt, _statins.txt
  python gwas_snp_annotator.py -i /path/to/snps/ -o annotated_results.csv \\
    --file-suffix _age.txt _sex.txt _statins.txt

  # Custom threshold and column names
  python gwas_snp_annotator.py -i /path/to/snps/ -p 1e-6 \\
    --snp-col RS_ID --chr-col CHROM --pos-col POS --p-col PVALUE

  # Custom flanking distance (50kb)
  python gwas_snp_annotator.py -i /path/to/snps/ -f 50000
        """
    )
    
    parser.add_argument('-i', '--input-dir', required=True, type=Path,
                       help='Directory containing SNP files (*.txt)')
    parser.add_argument('-o', '--output',
                       default='/Users/reneefonseca/Documents/UChicago/Dahl/Figures/Supp_tables/sig_loci/annotated_snps.csv',
                       help='Output CSV file path (default: Supp_tables/sig_loci/annotated_snps.csv)')
    parser.add_argument('-p', '--p-threshold', type=float, default=5e-8,
                       help='P-value threshold (default: 5e-8 Bonferroni)')
    parser.add_argument('-f', '--flank-distance', type=int, default=10000,
                       help='Flanking distance for gene search in bp (default: 10000)')
    
    # Column name arguments
    parser.add_argument('--snp-col', default='SNP',
                       help='SNP ID column name (default: SNP)')
    parser.add_argument('--chr-col', default='CHR',
                       help='Chromosome column name (default: CHR)')
    parser.add_argument('--pos-col', default='BP',
                       help='Position column name (default: BP)')
    parser.add_argument('--p-col', default='P',
                       help='P-value column name (default: P)')
    
    # File filtering
    parser.add_argument('--file-suffix', nargs='+', default=None,
                       help='Only process files ending with these suffixes (e.g., _age.txt _sex.txt _statins.txt)')
    
    args = parser.parse_args()
    
    # Validate input directory
    if not args.input_dir.exists():
        print(f"Error: Input directory does not exist: {args.input_dir}")
        sys.exit(1)
    
    # Initialize annotator
    annotator = SNPAnnotator(flank_distance=args.flank_distance)
    
    # Parse and filter SNPs
    print(f"Parsing SNP files from: {args.input_dir}")
    print(f"P-value threshold: {args.p_threshold}")
    print(f"Flanking distance: {args.flank_distance:,} bp")
    print("-" * 60)
    
    snp_df = annotator.parse_snp_files(
        args.input_dir, 
        args.p_threshold,
        snp_col=args.snp_col,
        chr_col=args.chr_col,
        pos_col=args.pos_col,
        p_col=args.p_col,
        file_patterns=args.file_suffix
    )
    
    # Annotate SNPs
    annotated_df = annotator.process_snps(
        snp_df,
        snp_col=args.snp_col,
        chr_col=args.chr_col,
        pos_col=args.pos_col,
        output_file=args.output
    )
    
    # Summary statistics
    print("\n" + "=" * 60)
    print("ANNOTATION SUMMARY")
    print("=" * 60)
    print(f"Total SNPs annotated: {len(annotated_df)}")
    print(f"SNPs with nearby genes (±{args.flank_distance:,}bp): {sum(annotated_df['nearest_gene'] != 'None within 10kb')}")
    print(f"SNPs within genes: {sum(annotated_df.get('gene_location', '') == 'within')}")
    print(f"SNPs in ClinVar: {sum(annotated_df.get('in_clinvar', False))}")
    
    if 'clinical_significance' in annotated_df.columns:
        clinvar_snps = annotated_df[annotated_df.get('in_clinvar', False)]
        if len(clinvar_snps) > 0:
            print("\nClinVar clinical significance breakdown:")
            sig_counts = clinvar_snps['clinical_significance'].value_counts()
            for sig, count in sig_counts.items():
                print(f"  {sig}: {count}")
    
    print("\n" + "=" * 60)


if __name__ == '__main__':
    main()
