rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "compiled_output.Rdata"))

setwd(paste0(fig_dir, "Supp_plots/PGSCvampPGS"))

pheno_list <- plot_phenos
alt_pgs    <- c("ampPGS", "PGSC_v_ampPGS")
clean_meth <- c("ampPGS", "PGSC vs ampPGS")
names(clean_meth) <- alt_pgs
base_cols  <- base_cols3[c('ampPGS', 'PGSC_v_pgs')]
names(base_cols) <- alt_pgs

# per-context summary for this PGSC-vs-ampPGS comparison
pgsc_amp_summary <- array(NA, dimnames = list(context_list, valid_pops, alt_pgs, c('r2','se')),
                          dim = c(length(context_list), length(valid_pops), length(alt_pgs), 2))

for(pop in valid_pops){
  for(context in context_list){

    output_loc <- output_valid[pheno_list,,, context, pop]

    # baseline R2
    base_r2  <- output_loc[, "ampPGS", "r2"]
    keep <- which( base_r2 > 0.01 )

    # % differences in R2
    update <- cbind(output_loc[,alt_pgs[1],"r2"],output_loc[,alt_pgs[2],"delta"])
    colnames(update) <- alt_pgs # the ampPGS col here is fake just so the plotting function works
    per_diff_r2 <- scale_by_pgs(update, base_r2)[keep,]
    SD          <- scale_by_pgs(output_loc[,alt_pgs,"sd"]   , base_r2)[keep,]
    CI25        <- scale_by_pgs(output_loc[,alt_pgs,"CI25"] , base_r2)[keep,]
    CI75        <- scale_by_pgs(output_loc[,alt_pgs,"CI75"] , base_r2)[keep,]

    avg_se      <- avg_se_fn(SD)
    m           <- build_plot_matrices(per_diff_r2, SD, CI25, CI75, sort_col = "PGSC_v_ampPGS", avg_se = avg_se)
    data_matrix <- m$data
    CI25_matrix <- m$CI25
    CI75_matrix <- m$CI75
    avg_r2s     <- data_matrix["avg_r2s", ]
    pgsc_amp_summary[context, pop, , ] <- cbind(avg_r2s, avg_se)
    # One-sample Wald z-test: H0 mean delta = 0, using same sandwich SE as the CIs
    pvals <- setNames(2 * pnorm(-abs(avg_r2s / avg_se)), alt_pgs)

    r2diffplot( data_matrix, CI25_matrix, CI75_matrix,
                filename=paste0("R2diff_PGSC_v_ampPGS_", context, '_', pop, ".png"),
                alt_pgs[2], clean_meth[2],
                pheno_cleaner(rownames(data_matrix)),
                context_mapping[[context]]$clean_context, pop, context )

    # print stuff:
    if(context == "sex" && pop == "white_euro"){
      print(paste0("P-values (avg R2% change vs 0) — ", context, " | ", pop))
      print(round(pvals, digits=2))
      print(pvals)
      print(paste0("# traits: ", length(data_matrix[-c(1,2),2])))
      print("avg R2 %Change: ")
      print(data_matrix["avg_r2s",])
      print("top PGSC R2 %Change: ")
      print(tail(data_matrix[,"PGSC_v_ampPGS"],n=4))
      print(paste0("# PGSC sig: ", length(CI25_matrix[CI25_matrix[,"PGSC_v_ampPGS"] > 0,"PGSC_v_ampPGS"])))
      print(paste0("# ampPGS sig: ", length(CI75_matrix[CI75_matrix[,"PGSC_v_ampPGS"] < 0,"PGSC_v_ampPGS"])))
      print("PGSC underperforms ampPGS: ")
      print(CI75_matrix[CI75_matrix[,"PGSC_v_ampPGS"] < 0,])

    }
  }
}

