rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "compiled_output.Rdata"))

setwd(paste0(fig_dir, "ED_Fig1"))

genie <- read.csv(paste0(home_dir,'GENIE/ali_ajhg_results.csv'), head=T )

match_list <- match_list[!names(match_list) %in% c(pheno_drop, pheno_sig_drop)]
match_list_cleaned <- setNames(match_list, sapply(names(match_list), pheno_cleaner))
pheno_match_vals <- unname(match_list_cleaned)

### for w1a w1b label axis as male v female or old v young
pheno_list <- plot_phenos
base_cols <- c(base_cols3, '#E88BB6')
names(base_cols) <- c(alt_pgs3,"1b")
w1a_labels <- c("males", "born 1952-1971", "statin users")
w1b_labels <- c("females", "born 1934-1951", "statin non-users")

taus <- array(NA,dim=c(length(pheno_list), length(context_list),2),
              dimnames = list(pheno_list, context_list, c('tau0','tau1')))
diff_r2 <- array(NA,dim=c(length(pheno_list), length(context_list), length(alt_pgs3)),
                     dimnames = list(pheno_list, context_list, alt_pgs3 ) )
per_diff_r2 <- array(NA,dim=c(length(pheno_list), length(context_list), length(alt_pgs3)),
                     dimnames = list(pheno_list, context_list, alt_pgs3 ) )
ws <- array(NA,dim=c(length(pheno_list), length(context_list),length(alt_pgs3),5),
            dimnames = list(pheno_list, context_list,alt_pgs3,c("w0","w0a","w0b","w1a","w1b")))
merged_list <- list()

for(context in context_list){
  for(pheno in pheno_list){
      diff_r2[pheno,context,alt_pgs3[2:3]] <- output_valid[pheno, alt_pgs3[2:3], 'delta',context, valid_pops[1]]
      per_diff_r2[pheno,context,alt_pgs3[2:3]] <- scale_by_pgs(output_valid[pheno, alt_pgs3[2:3], 'delta',context, valid_pops[1]],
                                                              output_valid[pheno, 'pgs', 'r2',context, valid_pops[1]])
      ws[pheno,context,alt_pgs3,"w0"] <- output_scaled[pheno, c('pgs','ampPGS','PGSC'),'w0', context]
      ws[pheno,context,"ampPGS",c('w0a','w0b')] <- output_scaled[pheno, 'ampPGS',c('w0a','w0b'), context]
      ws[pheno,context,"PGSC_v_pgs",c('w1a','w1b')] <- output_scaled[pheno, 'PGSC',c('w1a','w1b'), context]
      taus[pheno,context,'tau1'] <- output_valid[pheno, alt_pgs3[3], 'thresh', context, valid_pops[1]]
  }

  # set up genie & link to R2 diffs/ws
  g_con <- if (context == "statins") "statin" else context
  genie_context <- genie[genie[,'env'] == g_con,]
  R2s <- as.data.frame(cbind(diff_r2[pheno_list,context,alt_pgs3[2:3]],per_diff_r2[pheno_list,context,alt_pgs3[2:3]]))
  colnames(R2s) <- c(paste0("diff_",alt_pgs3[2:3]), alt_pgs3[2:3])
  PGSC_ws <- as.data.frame(ws[pheno_list,context,"PGSC_v_pgs",c("w0",'w1a','w1b')])
  R2s$clean_phenos <- pheno_cleaner(rownames(R2s))
  R2s$pheno_match <- pheno_match_vals
  PGSC_ws$pheno_match <- pheno_match_vals
  tau1_df <- data.frame(tau1 = taus[pheno_list, context, 'tau1'], pheno_match = pheno_match_vals)
  merged_df <- na.omit(merge(
    genie_context[, c('pheno','h2g','h2gxe','h2nxe')],
    cbind(R2s, PGSC_ws[, !names(PGSC_ws) %in% "pheno_match"],
          tau1_df[, !names(tau1_df) %in% "pheno_match"]),
    by.x = "pheno", by.y = "pheno_match", all.x = TRUE
  ))
  merged_list[[context]] <- merged_df
}

# build h2 gxe & het plot
panel_labels <- paste0("(", LETTERS[1:6], ")")

box_side   <- 6.8
mai_bottom <- 1.2; mai_left <- 1.5; mai_right <- 0.45
mai_top    <- 1.3
col_width  <- box_side + mai_left + mai_right
row_height <- box_side + mai_top + mai_bottom
title_width <- col_width - 0.4

cex_lab    <- 4
cex_axis   <- 3.5
cex_panel  <- 2
cex_cor    <- 1.9
cex_legend <- 3.5
cex_pts    <- 4
panel_lwd  <- 3.5

fit_cex_main <- function(label, target = title_width, cap = cex_lab, font = 1) {
  min(cap, target / strwidth(label, units = "inches", cex = 1, font = font))
}

draw_gxe_het_panel <- function(merged_df, xvals, xlab_str, panel_label,
                                context_label = NULL, show_ylab = FALSE,
                                show_legend = FALSE) {
  plot(xvals, merged_df[, "PGSC_v_pgs"],
       xlim = range(xvals), ylim = range(merged_df[, c("PGSC_v_pgs", "ampPGS")]),
       xlab = "", ylab = "", col = base_cols["PGSC_v_pgs"], pch = 16, cex = cex_pts, cex.axis = cex_axis)
  points(xvals, merged_df[, "ampPGS"], col = base_cols["ampPGS"], pch = 16, cex = cex_pts)
  abline(h = 0, col = "black")
  abline(v = 0, col = "black")
  title(xlab = xlab_str, cex.lab = cex_lab, line = 7)

  PGSC_cor   <- cor.test(merged_df[, "PGSC_v_pgs"], xvals, method = "pearson")
  ampPGS_cor <- cor.test(merged_df[, "ampPGS"],     xvals, method = "pearson")
  cor_line   <- paste(cor_label("PGSC", PGSC_cor), cor_label("ampPGS", ampPGS_cor), sep = "; ")

  if (!is.null(context_label)) {
    title(main = context_label, cex.main = fit_cex_main(context_label), font.main = 1, line = 5)
  }
  mtext(cor_line, side = 3, line = 1.2, cex = cex_cor, font = 1)

  if (show_ylab) {
    title(ylab = paste0("R\u00b2 %Change ", Ancs[1], " population"), cex.lab = cex_lab, line = 7)
    if (show_legend) {
      legend("right", legend = c("PGSC", "ampPGS"), bty = "n",
             col = base_cols[c("PGSC_v_pgs", "ampPGS")], pch = c(16, 16), cex = cex_legend)
    }
  }
}

png("R2perdiffvgxe_het.png", width = 3 * col_width, height = 2 * row_height, units = 'in', res = 300)
layout(matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE),
       widths = rep(col_width, 3), heights = c(row_height, row_height))

par(mai = c(mai_bottom, mai_left, mai_top, mai_right), lwd = panel_lwd, tcl=-0.8, mgp=c(3,2.5,0))
for (i in seq_along(context_list)) {
  context   <- context_list[i]
  merged_df <- merged_list[[context]]
  draw_gxe_het_panel(merged_df, merged_df$h2gxe, "GxC heritability", panel_labels[i],
                      context_label = context_mapping[[context]]$clean_context,
                      show_ylab = (context == "sex"), show_legend = TRUE)
}

par(mai = c(mai_bottom, mai_left, mai_top, mai_right), lwd = panel_lwd, tcl=-0.8, mgp=c(3,2.5,0))
for (i in seq_along(context_list)) {
  context   <- context_list[i]
  merged_df <- merged_list[[context]]
  draw_gxe_het_panel(merged_df, abs(merged_df$h2nxe), "Heteroskedasticity",
                      panel_labels[i + length(context_list)],
                      show_ylab = (context == "sex"))
}
dev.off()

png("Rhos_combined_we.png", width = 24, height = 16, units = 'in', res = 300)
par(mfrow = c(2, 3), mai = c(0.5,1,0.75,0.5))
# Top row: Rhos v gxe (GxC h2)
for (i in seq_along(context_list)) {
  context   <- context_list[i]
  merged_df <- merged_list[[context]]
  PGSC1a_rho <- merged_df$w1a/merged_df$w0
  PGSC1b_rho <- merged_df$w1b/merged_df$w0
  diff <- PGSC1a_rho - PGSC1b_rho
  plot(merged_df$h2gxe, diff,
       xlim = range(merged_df$h2gxe), ylim = range(diff),
       xlab = "", ylab = "",
       col = base_cols["PGSC_v_pgs"],
       pch=16, cex=4, cex.lab = 3, cex.axis=2.5)
  title(xlab = "GxC heritability", cex.lab = 3.5, line = 4)
  title(ylab = bquote("Diff(Gx" * .(str_to_sentence(context)) * ") " * rho[c]), cex.lab = 3.5, line = 4)
  abline(h=0, col="black")
  abline(v=0, col="black")
  mtext(panel_labels[i], side=3, adj=0.01, line = -2.8, cex=2, font=1)
  PGSC_cor <- cor.test(merged_df$h2gxe, diff, method = "pearson")
  legend('topright', bty = 'n', c(context_mapping[[context]]$clean_context,
                                   cor_label("", PGSC_cor)), cex = 3)
}
# Bottom row: Rhos v heteroskedasticity
for (i in seq_along(context_list)) {
  context   <- context_list[i]
  merged_df <- merged_list[[context]]
  PGSC1a_rho <- merged_df$w1a/merged_df$w0
  PGSC1b_rho <- merged_df$w1b/merged_df$w0
  avg <- rowMeans(cbind(PGSC1a_rho, PGSC1b_rho))
  plot(abs(merged_df$h2nxe), abs(avg),
       xlim = range(abs(merged_df$h2nxe)), ylim = range(abs(avg)),
       xlab = "", ylab = "",
       col = base_cols["PGSC_v_pgs"],
       pch=16, cex=4, cex.lab = 3, cex.axis=2.5)
  mtext(panel_labels[i + length(context_list)], side=3, adj=0.01, line = -2.8, cex=2, font=1)
  title(xlab = "Heteroskedasticity", cex.lab = 3.5, line = 4)
  title(ylab = bquote("Mean(Gx" * .(str_to_sentence(context)) * ") " * rho[c]), cex.lab = 3.5, line = 4)
  PGSC_cor <- cor.test(abs(merged_df$h2nxe), abs(avg), method = "pearson")
  legend('topright', bty = 'n', legend = cor_label("", PGSC_cor), cex = 3)
}
dev.off()





