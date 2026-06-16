# GxC SNP annotation

Annotates significant SNPs with the nearest gene (±10 kb, Ensembl GRCh37) and ClinVar clinical significance (MyVariant.info). Input SNP files come from `sig_loci_count.sh`.

## Setup

```bash
pip3 install pandas requests
```

## Run

```bash
# 1. Preview file format to confirm column names
python3 preview_snp_files.py /path/to/PGSC_data/sig_loci/best_snps/

# 2. Annotate (active contexts only)
python3 gwas_snp_annotator.py \
  -i /path/to/PGSC_data/sig_loci/best_snps/ \
  -o annotated_gwas_snps.csv \
  --snp-col ID --p-col P \
  --file-suffix _age.txt _sex.txt _statins.txt
```

If files contain only rsIDs (no CHR/BP), positions are looked up via Ensembl automatically. Runs at ~3–5 SNPs/sec (API rate limits); results are written incrementally, so a run can be resumed by re-running. Test on a subset first (`head -n 1001 file > test/file`). Run `--help` for all options.

## Options

| Flag | Default | Description |
|---|---|---|
| `-i, --input-dir` | required | Directory of SNP `.txt` files |
| `-o, --output` | `annotated_snps.csv` | Output CSV |
| `-p, --p-threshold` | `5e-8` | P-value threshold |
| `-f, --flank-distance` | `10000` | Gene search distance (bp) |
| `--snp-col` / `--p-col` | `SNP` / `P` | SNP-ID and p-value column names |
| `--chr-col` / `--pos-col` | `CHR` / `BP` | Optional position columns |
| `--file-suffix` | all `.txt` | Restrict to files ending with these suffixes |

## Output

The input columns plus: `nearest_gene`, `gene_distance`, `gene_location` (within/upstream/downstream), `gene_biotype`, `all_nearby_genes`, `in_clinvar`, `clinical_significance`, `conditions`, `clinvar_id`, `source_file`.
