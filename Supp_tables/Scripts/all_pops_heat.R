rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "summ_stats.Rdata"))
setwd(paste0(fig_dir, "Supp_tables"))

get_deltas <- function(phenos, pop, context)
  per_diff_r2_arr[phenos, 'PGSC_v_pgs', context, pop]

pair_cor <- function(phenos, pop1, pop2, context) {
  if (length(phenos) < 3) return(NA_real_)
  cor(get_deltas(phenos, pop1, context), get_deltas(phenos, pop2, context),
      method = "pearson", use = "pairwise.complete.obs")
}

pair_pval <- function(phenos, pop1, pop2, context) {
  if (length(phenos) < 3) return(NA_real_)
  cor.test(get_deltas(phenos, pop1, context), get_deltas(phenos, pop2, context),
           method = "pearson", exact = FALSE)$p.value
}

# Table: cor(WE deltas, pop deltas) per context x population, using the pairwise pheno set
alt_pops      <- valid_pops[c(2,3)]
alt_pop_names <- c("African", "Asian")

heat3_rows <- vector("list", length(context_list) * length(alt_pops))
idx <- 1
for(ctx_idx in seq_along(context_list)){
  context <- context_list[ctx_idx]
  pheno_we <- plot_phenos[which(pgs_r2_arr[plot_phenos, context, "white_euro"] > 0.01)]

  for(p_idx in seq_along(alt_pops)){
    pop      <- alt_pops[p_idx]
    pop_name <- alt_pop_names[p_idx]

    pheno_pop <- plot_phenos[which(pgs_r2_arr[plot_phenos, context, pop] > 0.01)]
    ph_pair   <- intersect(pheno_we, pheno_pop)

    heat3_rows[[idx]] <- data.frame(
      Context     = context,
      Population  = pop_name,
      Correlation = pair_cor(ph_pair,  "white_euro", pop, context),
      Pvalue      = pair_pval(ph_pair, "white_euro", pop, context),
      N           = length(ph_pair),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
}
heat3_data <- do.call(rbind, heat3_rows)

write.csv(heat3_data, "context_pop_correlation_table.csv", row.names = FALSE)


