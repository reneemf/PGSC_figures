rm(list = ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "compiled_output.Rdata"))

setwd(paste0(fig_dir, "Fig6"))
alt_pgs   <- alt_pgs2[1:2]          # ampPGS, PGSC_v_pgs
base_cols <- base_cols3[alt_pgs]

for(context in context_list) {
  base_r2  <- output_valid[scale_list, "pgs", "r2", context, "white_euro"]
  ampPGS_r2 <- as.data.frame(output_valid[scale_list, "ampPGS", "r2", context, "white_euro"])
  ampPGS_r2$clean_phenos <- pheno_cleaner(rownames(ampPGS_r2))
  ampPGS_r2$clean_phenos <- sub("^(.*) 0$", "log(\\1)", ampPGS_r2$clean_phenos)
  delta_r2 <- output_valid[scale_list, alt_pgs, "delta", context, "white_euro"]
  per_diff_r2 <- as.data.frame(scale_by_pgs(delta_r2, base_r2))
  per_diff_r2$clean_phenos <- pheno_cleaner(rownames(per_diff_r2))
  per_diff_r2$clean_phenos <- sub("^(.*) 0$", "log(\\1)", per_diff_r2$clean_phenos)

  SD <- as.data.frame(scale_by_pgs(output_valid[scale_list, alt_pgs, "sd", context, "white_euro"], base_r2))
  CI25 <- as.data.frame(scale_by_pgs(output_valid[scale_list, alt_pgs, "CI25", context, "white_euro"], base_r2))
  CI75 <- as.data.frame(scale_by_pgs(output_valid[scale_list, alt_pgs, "CI75", context, "white_euro"], base_r2))

  # Sort by phenotype order
  pheno_order <- c(
    "log(Testosterone)","Testosterone","log(SHBG)","SHBG",
    "log(Lipoprotein a)","Lipoprotein a","log(LDL)","LDL",
    "log(Height)","Height","log(HDL)","HDL",
    "log(HbA1c)","HbA1c",
    "log(Bilirubin)", "Bilirubin",
    "log(Arm fat-free mass left)","Arm fat-free mass left",
    "log(Apolipoprotein B)","Apolipoprotein B",
    "log(Apolipoprotein A)","Apolipoprotein A"
  )
  sort_index <- match(pheno_order, per_diff_r2$clean_phenos)
  valid_index <- which(!is.na(sort_index))

  ampPGS_r2 <- ampPGS_r2[sort_index[valid_index], , drop = FALSE]
  per_diff_r2 <- per_diff_r2[sort_index[valid_index], , drop = FALSE]
  SD <- SD[sort_index[valid_index], , drop = FALSE]
  CI25 <- CI25[sort_index[valid_index], , drop = FALSE]
  CI75 <- CI75[sort_index[valid_index], , drop = FALSE]

  # Split log vs non-log phenotypes
  is_log <- rownames(per_diff_r2) %in% pheno_drop

  #Summary stats with proper CIs
  r2_mean     <- colMeans(per_diff_r2[!is_log, c("ampPGS", "PGSC_v_pgs")]) 
  log_r2_mean <- colMeans(per_diff_r2[is_log, c("ampPGS", "PGSC_v_pgs")]) 

  subset_SD <- SD[!is_log, c("ampPGS", "PGSC_v_pgs")]
  log_subset_SD <- SD[is_log, c("ampPGS", "PGSC_v_pgs")]
  avg_se     <- avg_se_fn(subset_SD)
  avg_log_se <- avg_se_fn(log_subset_SD)

  # Build matrices for plotting
  data_matrix <- rbind(log_r2_mean, r2_mean, c(0, 0),
                       per_diff_r2[, c("ampPGS", "PGSC_v_pgs")])
  colnames(data_matrix) <- c("ampPGS", "PGSCgen")

  CI25_matrix <- rbind( log_r2_mean - 1.96*avg_log_se,
                        r2_mean     - 1.96*avg_se, c(0, 0),
                        CI25[, c("ampPGS", "PGSC_v_pgs")])

  CI75_matrix <- rbind(log_r2_mean + 1.96*avg_log_se,
                       r2_mean     + 1.96*avg_se, c(0, 0),
                       CI75[, c("ampPGS", "PGSC_v_pgs")])

  y_vals <- seq_len(nrow(data_matrix))
  plot_range <- range(data_matrix, na.rm = TRUE)
  bg_blocks <- ceiling(((y_vals - 3) / 2))
  bg_color <- c("white",detail_col["grey_1"])# "grey90")

  png(paste0("R2diff_scaled_",context,".png"), width = 10, height = 8, units = 'in', res = 300)
  par(mai = c(1, 2.5, 0.25, 0.25))
  plot(plot_range, range(y_vals), type = "n", xlab = "R\u00b2 %Change in European population", ylab = "",
       las = 1, xaxt = "s", yaxt = "n", cex.lab = 1.5, bty = "n")

  # Draw alternating color backgrounds
  for (i in seq_len(max(bg_blocks, na.rm = TRUE))) {
    idx <- which(bg_blocks == i)
    rect(par("usr")[1], min(idx) - 0.5, par("usr")[2], max(idx) + 0.5,
         col = bg_color[(i %% 2) + 1], border = NA)
  }

  # Points and CI arrows
  points(data_matrix[, "ampPGS"], y_vals, pch = 16,
         col = ifelse(data_matrix[, "ampPGS"] == 0, "white", base_cols["ampPGS"]), cex = 1.5)
  arrows(CI25_matrix[, "ampPGS"], y_vals, CI75_matrix[, "ampPGS"], y_vals,
         angle = 90, code = 3, lwd = 2, length = 0.05,
         col = ifelse(data_matrix[, "ampPGS"] == 0, "black", base_cols["ampPGS"]))

  points(data_matrix[, "PGSCgen"], y_vals, pch = 16,
         col = ifelse(data_matrix[, "PGSCgen"] == 0, "white", base_cols["PGSC_v_pgs"]), cex = 1.5)
  arrows(CI25_matrix[, "PGSC_v_pgs"], y_vals, CI75_matrix[, "PGSC_v_pgs"], y_vals,
         angle = 90, code = 3, lwd = 2, length = 0.05,
         col = ifelse(data_matrix[, "PGSCgen"] == 0, "black", base_cols["PGSC_v_pgs"]))

  # Axis and legends
  axis_labels <- c("Avg across log(phenos)", "Avg across phenos", " ", per_diff_r2$clean_phenos)
  axis(2, at = y_vals, labels = axis_labels, las = 1, cex.axis = 1)
  abline(v = 0, col = "black")
  abline(h = 3, col = "black")

  legend("topright", bty = "n", legend = context_mapping[[context]]$clean_context, cex = 1.5)
  if (context == "statins") {
    legend("bottomright", legend = c("PGSC", "ampPGS"), bty = "n",
           col = c(base_cols["PGSC_v_pgs"], base_cols["ampPGS"]), pch = 16, cex = 1.2, lwd = 2)
  }
  dev.off()
}