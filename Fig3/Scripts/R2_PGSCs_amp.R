rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "summ_stats.Rdata"))

setwd(paste0(fig_dir, "Fig3"))

base_cols <- base_cols3[c("ampPGS", "PGSC_v_pgs")]
port_avg_list <- c()

for (pop in valid_pops) {
  for (context in context_list) {

    keep <- which(pgs_r2_arr[, context, pop] > 0.01)

    per_diff_r2 <- per_diff_r2_arr[keep, , context, pop]
    SD          <- SD_pct_arr[keep,      , context, pop]
    CI25        <- CI25_pct_arr[keep,    , context, pop]
    CI75        <- CI75_pct_arr[keep,    , context, pop]

    avg_r2s <- context_summary[context, pop, , 'r2']
    avg_se  <- context_summary[context, pop, , 'se']
    pvals   <- pvals_vs0[context, pop, ]

    m           <- build_plot_matrices(per_diff_r2, SD, CI25, CI75,
                                       avg_r2s = avg_r2s, avg_se = avg_se)
    data_matrix <- m$data
    CI25_matrix <- m$CI25
    CI75_matrix <- m$CI75

    if(pop %in% valid_pops[2:3]){
      port_avg_list <- append(port_avg_list, data_matrix["avg_r2s", "PGSC_v_pgs"])
    }

    r2diffplot(data_matrix, CI25_matrix, CI75_matrix,
               filename  = paste0("R2diff_pgsc_amp_", context, '_', pop, ".png"),
               alt_pgs2[1:2], clean_meth2[1:2],
               pheno_cleaner(rownames(data_matrix)),
               context_mapping[[context]]$clean_context, pop, context)

  }
}
print(mean(port_avg_list))
