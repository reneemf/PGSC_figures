# PGSC figures

Scripts to reproduce the figures in the PGSC manuscript.

📄 **Manuscript:** [TITLE](MANUSCRIPT_URL) 

The pipeline that builds the underlying scores (GXC GWAS, PGSC, bootstrapped R², etc.) is in the companion method repository [here](repo_URL).

Three methods are compared throughout: **PGS** (additive baseline), **ampPGS** (genome-wide amplification), and **PGSC** (locus-specific GxC), across 48 UKB traits, three contexts (sex, age, statins), and three populations (European, African, Asian).

## Requirements

R ≥ 4.2:

```r
install.packages(c("ggplot2", "ggrepel", "gridExtra", "grid", "stringr", "plotrix", "dplyr"))
```

## Input data

UK Biobank and BioMe data are access-controlled, so this repo ships code plus the tables in `Supp_tables/` only. Supply your own outputs (produced by [PGSC_repo](repo_URL)) and update `data_dir` in [`config.R`](config.R) accordingly.

## Run

```bash
export PGSC_HOME="/path/to/Figures_repo/"   # repo root; figures write here
# edit data_dir in config.R, then:
Rscript process_data.R                       # writes Rdata/ + Supp_tables/summ_stats/*.csv
```
Then run any figure script (independent, any order).

## Figures

| Figure | Script(s) | 
|---|---|
| Fig2 — simulation panels | `Fig2/` (`run_sims.sh` → `collapse_results.R` → `panel_summary_plots.R`) | 
| Fig3 — per-phenotype R²% change | `Fig3/Scripts/R2_PGSCs_amp.R` | 
| Fig4 — cross-population portability | `Fig4/Scripts/all_pops_port.R` | 
| Fig5 — BioMe replication | `Fig5/Scripts/R2_PGSCs_amp_BioMe.R` | 
| Fig6 — log-transformed phenotypes | `Fig6/Scripts/R2_scaled_phenos_we.R` | 
| ED_Fig1 — GENIE GxC heritability comparison | `ED_Fig1/Scripts/R2s_ws_GENIE.R` | 
| ED_Fig2 — PGSC tuning parameters vs R²% change | `ED_Fig2/Scripts/ws_test_scaled.R` | 
| Supp_plots | `Supp_plots/Scripts/*.R` | 

Each script writes PNGs (and CSVs) into its own figure directory. Shared paths, phenotype lists, colors, and helper functions live in [`config.R`](config.R), sourced by every script.

## Supplementary tables

`Supp_tables/` holds `Phenotypes.xlsx` (static phenotype/context reference) plus scripts that generate tables:

| Script | Output |
|---|---|
| `process_data.R` (repo root) | `Supp_tables/summ_stats/summ_stats_{context}_{pop}.csv` | 
| `Supp_tables/Scripts/pop_counts.R` | `Supp_tables/pop_count/` count-table PNG | 
| `Supp_tables/Scripts/all_pops_corr.R` | `Supp_tables/context_pop_correlation_table.csv` 
| `Supp_tables/Scripts/sig_loci_analysis/gwas_snp_annotator.py` | `Supp_tables/sig_loci/annotated_gwas_snps.csv` | 

In the `summ_stats` CSVs, PGS columns are absolute R² while ampPGS/PGSC columns are ΔR² over baseline PGS; weight and context columns (`c_effect`, `pgs_gxc_cor`) come from the z-scored pipeline (context columns in `white_euro` files only).

Notes:
- **Fig2** simulations are self-contained and documented separately in [`Fig2/README.md`](Fig2/README.md).
- **ED_Fig1** compares PGSC results against the GENIE GxC heritability estimates from [Pazokitoroudi et al., *Am J Hum Genet* 2024](https://doi.org/10.1016/j.ajhg.2024.05.015).
- The GxC-significant **locus annotation** step has its own [README](Supp_tables/Scripts/sig_loci_analysis/README.md).

## License

MIT — see [LICENSE](LICENSE). Renée Fonseca, University of Chicago, Dahl Lab.
