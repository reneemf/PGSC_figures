rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "summ_stats.Rdata"))

setwd(paste0(fig_dir, "Supp_plots/C_effect_corr"))
input_dir <- paste0(home_dir, "r2_out/scaled/valid/")

# --- Define analysis parameters ---
types     <- alt_pgs2[1:2]          # ampPGS, PGSC_v_pgs
base_cols <- base_cols3[types]

# --- Build per_diff_r2 in long format ---
per_diff_r2_df <- do.call(rbind, lapply(context_list, function(context) {
  do.call(rbind, lapply(types, function(type) {
    data.frame(
      trait      = plot_phenos,
      context    = context,
      score_type = type,
      r2_diff    = as.numeric(per_diff_r2_arr[plot_phenos, type, context, "white_euro"]),
      stringsAsFactors = FALSE
    )
  }))
}))

per_diff_r2_df$clean_phenos <- pheno_cleaner(per_diff_r2_df$trait)

# --- Read p-values and effect sizes across all contexts ---
plot_data <- do.call(rbind, lapply(context_list, function(context) {
  pattern <- paste0("^misc_.*_white_euro_", context, "\\.txt$")
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  do.call(rbind, lapply(files, function(file) {
    df <- as.data.frame(t(read.table(file, header = TRUE, stringsAsFactors = FALSE)))
    trait <- gsub(paste0("misc_(.*)_white_euro_", context, "\\.txt"), "\\1", basename(file))
    df$trait <- trait
    df$context <- context
    df
  }))
}))

plot_data$clean_phenos <- pheno_cleaner(plot_data$trait)

# --- Read p-values and effect sizes across all contexts ---
merged_data <- merge(
  per_diff_r2_df,
  plot_data[, c("c_effect", "pgs_gxc_cor", "context", "clean_phenos")],
  by = c("context", "clean_phenos"),
  all.x = TRUE
)

# --- Split by score type ---
merged_PGSC <- merged_data[merged_data$score_type == "PGSC_v_pgs",]
merged_ampPGS <- merged_data[merged_data$score_type == "ampPGS",]

# Helper: drop testosterone for sex context, pass through otherwise
drop_testosterone <- function(df, con) {
  if (con == "sex") df[df$trait != "Testosterone674178", ] else df
}

# PGSC R2s vs effect of C
library(plotrix)

# break boundaries — adjust after checking:
# sort(merged_data[merged_data$context == "sex", "r2_diff"], decreasing = TRUE)[1:5]
break_lo <-35 #32.3 & amp 8
break_hi <- 100 #104.5 & amp: 105

# --- Helper: build gap.plot or standard plot depending on context ---
plot_panel <- function(x1, y1, x2, y2, ref_x, ref_y,
                       xlab, ylab, con, log_y = FALSE) {
  log_arg <- if (log_y) "y" else ""
  xlim    <- range(ref_x, na.rm = TRUE)
  ylim    <- range(ref_y, na.rm = TRUE)

  if (con == "sex") {
    gap.plot(x1, y1,
             gap = c(break_lo, break_hi), gap.axis = "x",
             brw = 0.025, breakcol="red", #white
             xlim = xlim, ylim = ylim,
             col = base_cols["PGSC_v_pgs"], pch = 16,
             xlab = xlab, ylab = ylab, #log = log_arg, # not working in gap.plot
             main = '')
    points(x2, y2, col = base_cols["ampPGS"], pch = 16)
  } else {
    plot(x1, y1,
         xlim = xlim, ylim = ylim,
         col = base_cols["PGSC_v_pgs"], pch = 16,
         xlab = xlab, ylab = ylab,
         log = log_arg, main = '')
    points(x2, y2, col = base_cols["ampPGS"], pch = 16)
  }
}


# --- PNG output ---
png("R2_vs_corr_Ceffect.png", width = 9, height = 6, units = 'in', res = 300)
par(mfrow = c(2, 3), mar = c(4, 4, 0.5, 0.5) + 0.1)
panel_labels <- paste0("(",LETTERS[1:6],")")
# Row 1: abs(main effect of C)
for (i in seq_along(context_list)) {
  con       <- context_list[i]
  pgsc_con  <- merged_PGSC[merged_PGSC$context == con, ]
  ampgs_con <- merged_ampPGS[merged_ampPGS$context == con, ]
  ref_con   <- merged_data[merged_data$context == con, ]

  plot_panel(pgsc_con$r2_diff,  abs(pgsc_con$c_effect),
             ampgs_con$r2_diff, abs(ampgs_con$c_effect),
             ref_con$r2_diff,   abs(ref_con$c_effect),
             xlab = "", ylab = "",
             con = con)
  abline(v = 0)
  mtext(panel_labels[i], side=3, adj=0.008, line = -1, cex=0.7, font=1)
  title(xlab = "R\u00b2 % change", cex.lab = 1.2, line = 2.5)
  if(con == "sex"){title(ylab = "abs(main effect of C)", cex.lab = 1.2, line = 2.5)}
  PGSC_cor   <- cor.test(pgsc_con$r2_diff,  abs(pgsc_con$c_effect),  method = "pearson")
  ampPGS_cor <- cor.test(ampgs_con$r2_diff, abs(ampgs_con$c_effect), method = "pearson")
  legend('topright', bty = 'n',
         c(context_mapping[[con]]$clean_context,
           cor_label("PGSC",   PGSC_cor),
           cor_label("ampPGS", ampPGS_cor)),
         cex = 1)
  if (con == "statins") {
    legend("right", legend = c("PGSC", "ampPGS"), bty = "n",
           col = c(base_cols["PGSC_v_pgs"], base_cols["ampPGS"]), pch = 16, cex = 1)
  }
}

# Row 2: corr(pgs,gxc)^2
for (i in seq_along(context_list)) {
  con       <- context_list[i]
  pgsc_con  <- merged_PGSC[merged_PGSC$context == con, ]
  ampgs_con <- merged_ampPGS[merged_ampPGS$context == con, ]
  ref_con   <- merged_data[merged_data$context == con, ]

  plot_panel(pgsc_con$r2_diff,       pgsc_con$pgs_gxc_cor^2,
             ampgs_con$r2_diff,      ampgs_con$pgs_gxc_cor^2,
             ref_con$r2_diff,        ref_con$pgs_gxc_cor^2,
             xlab = "", ylab = "",
             con = con)
  abline(v = 0)
  mtext(panel_labels[i + length(context_list)], side=3, adj=0.008, line = -1, cex=0.7, font=1)
  title(xlab = "R\u00b2 % change", cex.lab = 1.2, line = 2.5)
  if(con == "sex"){title(ylab = expression("corr(pgs,gxc)"^2), cex.lab = 1.2, line = 2.5)}
  PGSC_cor   <- cor.test(pgsc_con$r2_diff,  pgsc_con$pgs_gxc_cor^2,  method = "pearson")
  ampPGS_cor <- cor.test(ampgs_con$r2_diff, ampgs_con$pgs_gxc_cor^2, method = "pearson")
  legend('topright', bty = 'n',
         c(context_mapping[[con]]$clean_context,
           cor_label("PGSC",   PGSC_cor),
           cor_label("ampPGS", ampPGS_cor)),
         cex = 1)
}

dev.off()

