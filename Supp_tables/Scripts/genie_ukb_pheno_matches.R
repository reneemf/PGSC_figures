rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))

# ── GENIE ↔ UKB phenotype match table ─────────────────────────────────────────
# Mirrors the phenotype selection in ED_Fig1/Scripts/R2s_ws_GENIE.R: the GENIE
# phenotypes actually used in analysis are the UKB phenotypes that (1) survive the
# pheno_drop / pheno_sig_drop filtering, (2) have a non-NA GENIE match in
# match_list, and (3) appear in the GENIE results file. Writes the pairing to a
# CSV in Supp_tables/.

load(paste0(rdata_dir, "compiled_output.Rdata"))   # provides plot_phenos, pheno_sig_drop

genie <- read.csv(paste0(home_dir, "GENIE/ali_ajhg_results.csv"), header = TRUE)

# drop phenotypes not carried into the analysis (same filter as R2s_ws_GENIE.R)
match_list_used <- match_list[!names(match_list) %in% c(pheno_drop, pheno_sig_drop)]

# keep only UKB phenotypes with a GENIE match that is present in the GENIE file
match_list_used <- match_list_used[!is.na(match_list_used) &
                                     match_list_used %in% unique(genie$pheno)]

pheno_table <- data.frame(
  PGSC_phenotype       = pheno_cleaner(names(match_list_used)),
  GENIE_phenotype     = unname(match_list_used),
  row.names           = NULL,
  stringsAsFactors    = FALSE
)
pheno_table <- pheno_table[order(pheno_table$PGSC_phenotype), ]

out_file <- paste0(fig_dir, "Supp_tables/genie_pgsc_pheno_matches.csv")
write.csv(pheno_table, file = out_file, row.names = FALSE)

