rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "summ_stats.Rdata"))

setwd(paste0(fig_dir, "Supp_plots/loci_count"))

alt_pgs   <- type_list_valid[3:4]   # ampPGS, PGSC_v_pgs
base_cols <- base_cols3[alt_pgs]

for(context in context_list){
  clean_context <- context_mapping[[context]]$clean_context
  sig_file  <- paste0(home_dir, "sig_loci/sig_loci_", context, ".txt")
  sig_loci  <- read.table(sig_file)
  row.names(sig_loci) <- sig_loci[,1]
  per_diff_r2 <- na.omit(as.data.frame(per_diff_r2_arr[plot_phenos, alt_pgs2, context, "white_euro"]))
  per_diff_r2$clean_phenos <- pheno_cleaner(rownames(per_diff_r2))

  #subset datasets
  pheno_subset <- intersect(sig_loci[,1], row.names(per_diff_r2))
  sig_loci <- sig_loci[pheno_subset,]
  per_diff_r2 <- per_diff_r2[pheno_subset,]

  sig_log <- sig_loci[,2] + 1

  # R2 vs log(loci_count)
  png(paste0("R2vloci_countlog2_",context,".png"),
      width = 8, height = 8, units = 'in', res = 300)
  par(mai=c(1,1,0.5,0.5))
  # plot PGSC points
  plot(sig_log, per_diff_r2[,"PGSC_v_pgs"], log = 'x',
       xlim=range(sig_log,na.rm =T) + c(0,0.1) * diff(range(sig_log, na.rm =T)),
       ylim=range(per_diff_r2[,c(1:2)], na.rm =T) + c(-0.02,0.07) * diff(range(per_diff_r2[,c(1:2)], na.rm =T)),
       xlab="log(# sig GxC GWAS loci)", ylab="R\u00b2 %Change",
       col=base_cols["PGSC_v_pgs"], pch=16,cex=1.5,cex.lab = 1.5)
  # add ampPGS points
  points(sig_log,per_diff_r2[,"ampPGS"], col=base_cols["ampPGS"], pch=16,cex=1.5)
  # annotate top 2 PGSC values
  pgsc_idx <- order(per_diff_r2[,"PGSC_v_pgs"], decreasing=TRUE)[1:2]
  pgsc_names <- pheno_cleaner(rownames(per_diff_r2)[pgsc_idx])
  text(sig_log[pgsc_idx], per_diff_r2[pgsc_idx,"PGSC_v_pgs"],
       labels=pgsc_names, pos=1, offset=0.7, cex=1, col=base_cols["PGSC_v_pgs"])
  # annotate top 2 ampPGS values
  amp_idx <- order(per_diff_r2[,"ampPGS"], decreasing=TRUE)[1:2]
  amp_names <- pheno_cleaner(rownames(per_diff_r2)[amp_idx])
  text(sig_log[amp_idx], per_diff_r2[amp_idx,"ampPGS"],
       labels=amp_names, adj = c(0.4,-1.45), cex=1, col=base_cols["ampPGS"])
  # legends
  if(context == "sex"){
    legend("bottomright", legend = c("PGSC","amplification PGS"), bty = "n",
           col = c(base_cols["PGSC_v_pgs"], base_cols["ampPGS"]), pch = c(16, 16), cex = 1.3)
  }
  PGSC_cor   <- cor.test(per_diff_r2[,"PGSC_v_pgs"],  sig_log,  method = "pearson")
  ampPGS_cor <- cor.test(per_diff_r2[,"ampPGS"], sig_log, method = "pearson")
  legend('right', bty = 'n',
         c(context_mapping[[context]]$clean_context,
           cor_label("PGSC",   PGSC_cor),
           cor_label("ampPGS", ampPGS_cor)),
         cex = 1)

  dev.off()

}


