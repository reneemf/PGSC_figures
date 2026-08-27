rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "summ_stats.Rdata"))

setwd(paste0(fig_dir, "Supp_plots/C_effect_corr"))
input_dir <- paste0(home_dir, "r2_out/scaled/valid/")

types     <- alt_pgs2[1:2]          # ampPGS, PGSC_v_pgs
base_cols <- base_cols3[types]

# build per_diff_r2 in long format
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

# read p-values and effect sizes across all contexts
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

merged_data <- merge(
  per_diff_r2_df,
  plot_data[, c("c_effect", "pgs_gxc_cor", "context", "clean_phenos")],
  by = c("context", "clean_phenos"),
  all.x = TRUE
)

# split by score type
merged_PGSC <- merged_data[merged_data$score_type == "PGSC_v_pgs",]
merged_ampPGS <- merged_data[merged_data$score_type == "ampPGS",]

library(plotrix)

# break boundaries
break_lo <-35 
break_hi <- 100 

# build gap.plot or standard plot depending on context
plot_panel <- function(x1, y1, x2, y2, ref_x, ref_y,
                       xlab, ylab, con, log_y = FALSE) {
  log_arg <- if (log_y) "y" else ""
  xlim    <- range(ref_x, na.rm = TRUE)
  ylim    <- range(ref_y, na.rm = TRUE)

  if (con == "sex") {
    gap.plot(x1, y1,
             gap = c(break_lo, break_hi), gap.axis = "x",
             brw = 0.025, breakcol="red", 
             xlim = xlim, ylim = ylim,
             col = base_cols["PGSC_v_pgs"], pch = 16,
             xlab = xlab, ylab = ylab, 
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


# margin set up
measure_mai <- function(mar) {
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = 9, height = 6, units = 'in', res = 100)
  par(mar = mar)
  mai <- par("mai")
  dev.off()
  unlink(tmp)
  mai
}

# panel set up
panel_mai     <- measure_mai(c(4, 4, 3.2, 0.5) + 0.1)  
panel_size_in <- 3 - panel_mai[1] - panel_mai[3]        
col_width_in  <- panel_mai[2] + panel_size_in + panel_mai[4]
gap_in        <- 0.5 / 2.54 # 0.5cm column gap

row1_mar <- c(3, 3.5, 2.2, 0) + 0.1 # bottom, left, top, right (lines)
row2_mar <- c(3, 3.5, 2.2, 0) + 0.1 # tighter top margin (cor/pv line only)
row1_mai <- measure_mai(row1_mar)
row2_mai <- measure_mai(row2_mar)
row1_height_in <- row1_mai[1] + panel_size_in + row1_mai[3]
row2_height_in <- row2_mai[1] + panel_size_in + row2_mai[3]

fig_width_in  <- 3 * col_width_in + 2 * gap_in
fig_height_in <- row1_height_in + row2_height_in

# draw one row of 3 panels
panel_labels <- paste0("(", LETTERS[1:6], ")")
draw_row <- function(row_offset, mar, y_of, ylab_text, show_titles) {
  par(mar = mar, pty = "s")
  for (i in seq_along(context_list)) {
    con       <- context_list[i]
    pgsc_con  <- merged_PGSC[merged_PGSC$context == con, ]
    ampgs_con <- merged_ampPGS[merged_ampPGS$context == con, ]
    ref_con   <- merged_data[merged_data$context == con, ]

    plot_panel(pgsc_con$r2_diff,  y_of(pgsc_con),
               ampgs_con$r2_diff, y_of(ampgs_con),
               ref_con$r2_diff,   y_of(ref_con),
               xlab = "", ylab = "", con = con)
    abline(v = 0)
    if (show_titles) {
      title(main = context_mapping[[con]]$clean_context, font.main = 1, cex.main = 1.1, line = 1.5)
    }
    title(xlab = "R\u00b2 % change", cex.lab = 1.2, line = 2)
    if (con == "sex") title(ylab = ylab_text, cex.lab = 1.2, line = 2.5)

    PGSC_cor   <- cor.test(pgsc_con$r2_diff,  y_of(pgsc_con),  method = "pearson")
    ampPGS_cor <- cor.test(ampgs_con$r2_diff, y_of(ampgs_con), method = "pearson")
    mtext(paste(cor_label("PGSC", PGSC_cor), cor_label("ampPGS", ampPGS_cor), sep = "; "),
          side = 3, line = 0.3, cex = 0.5)

    if (show_titles && con == "statins") {
      legend("right", legend = c("PGSC", "ampPGS"), bty = "n",
             col = c(base_cols["PGSC_v_pgs"], base_cols["ampPGS"]), pch = 16, cex = 1)
    }
  }
}

png("R2_vs_corr_Ceffect.png", width = fig_width_in, height = fig_height_in, units = 'in', res = 300)
layout(matrix(c(1, 0, 2, 0, 3,
                4, 0, 5, 0, 6), nrow = 2, byrow = TRUE),
       widths  = c(col_width_in, gap_in, col_width_in, gap_in, col_width_in),
       heights = c(row1_height_in, row2_height_in))

draw_row(0, row1_mar, function(df) abs(df$c_effect), "abs(main effect of C)", show_titles = TRUE)
draw_row(3, row2_mar, function(df) df$pgs_gxc_cor^2, expression("corr(pgs,gxc)"^2), show_titles = FALSE)

dev.off()

