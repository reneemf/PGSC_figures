rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "compiled_output.Rdata"))

setwd(paste0(fig_dir, "Supp_plots/pgs"))

for(pop in valid_pops){
  files <- list.files(paste0(home_dir, "r2_out/valid"), pattern = paste0("^misc_.*\\_", pop, "_sex.txt$"), full.names = TRUE)
  plot_data <- data.frame()
  for (file in files) {
    df <- as.data.frame(t(read.table(file, header = TRUE, stringsAsFactors = FALSE)))
    trait <- gsub(paste0("misc_(.*)_", pop, "_sex.txt"), "\\1", basename(file))
    df$trait <- trait
    plot_data <- rbind(plot_data, df)
  }
  plot_data$clean_phenos <- pheno_cleaner(plot_data$trait)
  save(plot_data, file = paste0(rdata_dir, "misc_GxSex_pgs_output_", pop, ".Rdata"))

  # set up data sets
  pheno_list_update <- setdiff(pheno_list,pheno_drop)
  pgs_r2 <- output_valid[pheno_list_update,'pgs','r2',"sex",pop]
  clean_context <- context_mapping[["sex"]]$clean_context
  clean_phenos <- pheno_cleaner(names(pgs_r2))
  pgs_r2 <- as.data.frame(cbind(pgs_r2,clean_phenos))
  pgs_r2_pv <- merge(pgs_r2,plot_data)
  pv_sort <- na.omit(pgs_r2_pv[order(as.numeric(pgs_r2_pv$pgs_pv),decreasing=T),])
  r2_sort <- na.omit(pgs_r2_pv[order(as.numeric(pgs_r2_pv$pgs_r2),decreasing=F),])

  # r2 plot w pgs thresh
  r2_vals <- as.numeric(r2_sort$pgs_r2)
  png(paste0("pgsR2_pvthresh_",pop,".png"), width = 10, height = 15, units = 'in', res = 300)
  par(mai = c(1, 2.95, 0.25, 0.25))
  plot(x = r2_vals, y = seq_along(r2_vals), type = "p",
       pch = ifelse(r2_vals < 0.01, 21, 16),
       col = base_cols3["pgs"],
       xlab = paste0("PGS R\u00b2 in ", Ancs[pop], " population"), ylab = "",
       yaxt = "n", las = 1, cex.lab = 2, cex = 1.5, cex.axis = 1.5,
       xlim = range(r2_vals))
  axis(2, at = seq_along(r2_vals), labels = r2_sort$clean_phenos,
       las = 1, cex.axis = 1.3)
  abline(v = 0.01, col = "black")
  legend('bottomright', bty = 'n', clean_context, cex = 1.5)
  if(pop == "white_euro"){
    legend("right", legend = c("R\u00b2 < 1%", "R\u00b2 > 1%"), bty = 'n',
           col = base_cols3["pgs"], pch = c(21, 16), cex = 1.5)
  }
  dev.off()

  # drop list:
  print(r2_sort[r2_vals < 0.01, "clean_phenos"])

}
