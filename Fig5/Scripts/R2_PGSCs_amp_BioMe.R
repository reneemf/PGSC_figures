rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "summ_stats.Rdata"))

setwd(paste0(fig_dir, "Fig5"))

alt_pgs       <- alt_pgs3                   # c("pgs", "ampPGS", "PGSC_v_pgs")
base_cols     <- base_cols3
context       <- context_list[1]            # just working with sex here
clean_context <- context_mapping[[context]]$clean_context

# plot pgs R2s (all phenotypes, grey if R2 < 0.01)
pgs_r2_sort <- sort(biome_pgs_r2_all, na.last = NA)
png("R2_pgs_BioMe.png", width = 10, height = 12, units = 'in', res = 300)
par(mai = c(1, 2.95, 0.25, 0.25))
plot(x = pgs_r2_sort, y = seq_along(pgs_r2_sort),
     type = "p", pch = ifelse(pgs_r2_sort < 0.01, 21, 16),
     col = base_cols["pgs"],
     xlab = "PGS R\u00b2 in MSM population", ylab = "",
     yaxt = "n", las = 1, cex.lab = 2, cex = 1.5, cex.axis = 1.5,
     xlim = range(pgs_r2_sort))
axis(2, at = seq_along(pgs_r2_sort),
     labels = pheno_cleaner(names(pgs_r2_sort)),
     las = 1, cex.axis = 1.3)
abline(v = 0.01, col = "black")
legend('bottomright', bty = 'n', clean_context, cex = 1.5)
legend("right", legend = c("R\u00b2 > 1%", "R\u00b2 < 1%"), bty = 'n',
       col = base_cols["pgs"], pch = c(16, 21), cex = 1.5)
dev.off()

# Loop over all and cropped phenotype sets
pheno_sets <- c("all", "crop")
for (ver in pheno_sets) {
  if (ver == "all") {
    per_diff_r2 <- biome_per_diff_all
    SD          <- biome_sd_pct_all
    CI25        <- biome_CI25_all
    CI75        <- biome_CI75_all
  } else {
    per_diff_r2 <- biome_per_diff_crop
    SD          <- biome_sd_pct_crop
    CI25        <- biome_CI25_crop
    CI75        <- biome_CI75_crop
  }

  zero_set    <- rownames(per_diff_r2)[rowSums(per_diff_r2 == 0, na.rm = TRUE) > 0]
  avg_se      <- avg_se_fn(SD)
  m           <- build_plot_matrices(per_diff_r2, SD, CI25, CI75, avg_se = avg_se)
  data_matrix <- m$data
  data_matrix[zero_set, "PGSC_v_pgs"] <- 1E-10  # fix plotting issue where 0 is a white dot
  CI25_matrix <- m$CI25
  CI75_matrix <- m$CI75
  avg_r2s     <- data_matrix["avg_r2s", ]

  # One-sample Wald z-test: H0 mean delta = 0, using same sandwich SE as the CIs
  pvals <- setNames(2 * pnorm(-abs(avg_r2s / avg_se)), alt_pgs[2:3])

  # plot
  r2diffplot(data_matrix, CI25_matrix, CI75_matrix,
             filename  = paste0("R2diff_pgsc_amp_BioMe_", ver, ".png"),
             alt_pgs[2:3], clean_meth3[2:3],
             (pheno_cleaner(rownames(data_matrix))),
             clean_context, "BioMe", context)

  if (ver == "crop"){
    print("P-values (avg R2% change vs 0)")
    print(pvals)
    print("avg R2 %Change: ")
    print(data_matrix["avg_r2s", ])
    print(paste0("# phenos: ", length(data_matrix[-c(1,2),"PGSC_v_pgs"])))
    print("top PGSC R2 %Change: ")
    print(tail(data_matrix[, "PGSC_v_pgs"], n = 4))
    print("PGSC winners: ")
    print(CI25_matrix[CI25_matrix[, "PGSC_v_pgs"] > 0, ])
    print("ampPGS winners: ")
    print(CI25_matrix[CI25_matrix[, "ampPGS"] > 0, ])
  }
}
