rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "summ_stats.Rdata"))

setwd(paste0(fig_dir, "Supp_plots/loci_count"))

alt_pgs   <- alt_pgs2   # ampPGS, PGSC_v_pgs
base_cols <- base_cols3[alt_pgs]

# Annotate one point: thin leader line from the point to its label, placed a
# fixed vertical offset away (dir = -1 below, +1 above). An edge point's label
# is horizontally justified inward (right-most -> extends left, left-most ->
# extends right) so it can't run off the axis; interior points stay centered.
annotate_pt <- function(x, y, lab, col, dir, xall, y_span) {
  lab_y <- y + dir * 0.05 * y_span
  segments(x, y, x, lab_y, col = "black", lwd = 0.75)
  hj <- if (x == max(xall, na.rm = TRUE)) 1 else if (x == min(xall, na.rm = TRUE)) -0.01 else 0.5
  vj <- if (dir < 0) 1 else 0   # below the point -> top-aligned; above -> bottom-aligned
  text(x, lab_y, labels = lab, adj = c(hj, vj), cex = 1.2, col = col)
}

# 1x3 panel figure: one row, the three contexts (gxsex, gxage, gxstatins)
png("R2vloci_countlog2_1x3.png",
    width = 24, height = 8, units = 'in', res = 300)
par(mfrow = c(1, 3))

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

  # top margin holds the two title lines (context; cors/pvs). cex=1 undoes the
  # automatic mfrow text shrink so each 8x8 panel matches the old single plots.
  # only the gxsex panel gets the y-axis title.
  par(mai = c(1, 1, 1.35, 0.1), cex = 1, cex.axis = 1.5)
  ylab_use <- if(context == "sex") "R\u00b2 %Change" else ""

  # plot PGSC points
  plot(sig_log, per_diff_r2[,"PGSC_v_pgs"], log = 'x',
       xlim=range(sig_log,na.rm =T) + c(0,0.1) * diff(range(sig_log, na.rm =T)),
       ylim=range(per_diff_r2[,c(1:2)], na.rm =T) + c(-0.02,0.07) * diff(range(per_diff_r2[,c(1:2)], na.rm =T)),
       xlab="log(# sig GxC GWAS loci)", ylab=ylab_use,
       col=base_cols["PGSC_v_pgs"], pch=16,cex=1.5,cex.lab = 1.9)
  # add ampPGS points
  points(sig_log,per_diff_r2[,"ampPGS"], col=base_cols["ampPGS"], pch=16,cex=1.5)
  y_span <- diff(range(per_diff_r2[,c(1:2)], na.rm=TRUE))   # vertical offset scale for labels/leaders
  # annotate top 2 PGSC values below their points, top 2 ampPGS above theirs;
  # each label is tied to its point by a thin leader and justified inward at the axis edges
  pgsc_idx <- order(per_diff_r2[,"PGSC_v_pgs"], decreasing=TRUE)[1:2]
  pgsc_names <- pheno_cleaner(rownames(per_diff_r2)[pgsc_idx])
  for (k in seq_along(pgsc_idx))
    annotate_pt(sig_log[pgsc_idx[k]], per_diff_r2[pgsc_idx[k],"PGSC_v_pgs"],
                pgsc_names[k], base_cols["PGSC_v_pgs"], dir=-1, xall=sig_log, y_span=y_span)
  amp_idx <- order(per_diff_r2[,"ampPGS"], decreasing=TRUE)[1:2]
  amp_names <- pheno_cleaner(rownames(per_diff_r2)[amp_idx])
  for (k in seq_along(amp_idx))
    annotate_pt(sig_log[amp_idx[k]], per_diff_r2[amp_idx[k],"ampPGS"],
                amp_names[k], base_cols["ampPGS"], dir=1, xall=sig_log, y_span=y_span)

  # method legend: center-right of the gxsex panel only
  if(context == "sex"){
    legend("right", legend = c("PGSC","ampPGS"), bty = "n",
           col = c(base_cols["PGSC_v_pgs"], base_cols["ampPGS"]), pch = c(16, 16), cex = 1.3)
  }

  # context + correlations in the top title region:
  #   line 1 = context, line 2 = PGSC (left) and ampPGS (right) cor & pv
  PGSC_cor   <- cor.test(per_diff_r2[,"PGSC_v_pgs"], sig_log, method = "pearson")
  ampPGS_cor <- cor.test(per_diff_r2[,"ampPGS"],     sig_log, method = "pearson")

  # report the Sex-context correlations between log(# sig GxCWAS loci) and R2 %Change
  if(context == "sex"){
    cat("\n=== Context: Sex \u2014 cor(log(# sig GxCWAS loci), R2 %Change) ===\n")
    cat(sprintf("PGSC:   r = %.4f, p = %.4g\n", PGSC_cor$estimate,   PGSC_cor$p.value))
    cat(sprintf("ampPGS: r = %.4f, p = %.4g\n", ampPGS_cor$estimate, ampPGS_cor$p.value))

    # same correlations, excluding testosterone from the calculation
    keep <- rownames(per_diff_r2) != "Testosterone674178"
    PGSC_cor_noT   <- cor.test(per_diff_r2[keep,"PGSC_v_pgs"], sig_log[keep], method = "pearson")
    ampPGS_cor_noT <- cor.test(per_diff_r2[keep,"ampPGS"],     sig_log[keep], method = "pearson")
    cat("--- excluding testosterone ---\n")
    cat(sprintf("PGSC:   r = %.4f, p = %.4g\n", PGSC_cor_noT$estimate,   PGSC_cor_noT$p.value))
    cat(sprintf("ampPGS: r = %.4f, p = %.4g\n", ampPGS_cor_noT$estimate, ampPGS_cor_noT$p.value))
  }

  mtext(clean_context, side = 3, line = 3.85, cex = 1.8)
  if(context == "sex"){
    # gxsex panel: two cor/pv lines -- full analysis, then excluding testosterone.
    # the excluding-testosterone p is rounded to 3 dp so the small ampPGS p survives.
    mtext(cor_label("PGSC",   PGSC_cor),   side = 3, line = 2.3, adj = 0, cex = 1.4)
    mtext(cor_label("ampPGS", ampPGS_cor), side = 3, line = 2.3, adj = 1, cex = 1.4)

    # "excluding testosterone" caption flanked by drawn rules (in place of --- dashes).
    # y is a margin line converted to user units so the rules sit level with the text.
    # x-axis is log10, so par("usr") x is in log10 units: work there, then 10^ back
    # to data scale for text()/segments().
    usr <- par("usr")
    y_cap   <- usr[4] + 1.75 * par("csi") * (usr[4] - usr[3]) / par("pin")[2]
    cap     <- "excluding testosterone"
    xc      <- mean(usr[1:2])
    half    <- strwidth(cap, cex = 1.2) / 2
    gap     <- strwidth("n", cex = 1.2)
    seg_len <- 0.10 * (usr[2] - usr[1])
    text(10^xc, y_cap, cap, cex = 1.2, xpd = NA)
    segments(10^(xc - half - gap - seg_len), y_cap, 10^(xc - half - gap), y_cap, lwd = 1.5, xpd = NA)
    segments(10^(xc + half + gap), y_cap, 10^(xc + half + gap + seg_len), y_cap, lwd = 1.5, xpd = NA)

    noT_label <- function(name, cr)
      paste0(name, " Cor: ", round(cr$estimate, 2), ", p: ", round(cr$p.value, 3))
    mtext(noT_label("PGSC",   PGSC_cor_noT),   side = 3, line = 0.2, adj = 0, cex = 1.4)
    mtext(noT_label("ampPGS", ampPGS_cor_noT), side = 3, line = 0.2, adj = 1, cex = 1.4)
  } else {
    mtext(cor_label("PGSC",   PGSC_cor),   side = 3, line = 2.3, adj = 0, cex = 1.4)
    mtext(cor_label("ampPGS", ampPGS_cor), side = 3, line = 2.3, adj = 1, cex = 1.4)
  }
}

dev.off()


