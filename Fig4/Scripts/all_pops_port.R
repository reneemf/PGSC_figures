rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "summ_stats.Rdata"))
library(ggplot2)
library(ggrepel)

setwd(paste0(fig_dir, "Fig4"))

port_comps <- list(
  list(i = 1, j = 2, pv = "pgs_amp",  lvl = 1),
  list(i = 2, j = 3, pv = "amp_pgsc", lvl = 2),
  list(i = 1, j = 3, pv = "pgs_pgsc", lvl = 3)
)
tick_len <- 0.015
text_off <- 0.01
brk_step <- 0.05

png("avg_portability.png", width = 10, height = 6, units = "in", res = 300)
par(mfrow = c(1,2), mar = c(4,5,2,1))
for (anc in valid_pops[2:3]) {
  means    <- t(avg_port[, anc, , 'mean'])
  upper    <- t(avg_port[, anc, , 'mean'] + 1.96 * avg_port[, anc, , 'se'])
  lower    <- t(avg_port[, anc, , 'mean'] - 1.96 * avg_port[, anc, , 'se'])
  ylim_top <- max(upper, na.rm=TRUE) + 3 * brk_step + text_off + 0.03
  bp <- barplot(means, beside = TRUE, col = base_cols3[c('pgs','ampPGS','PGSC_v_pgs')],
                las = 1, ylim = c(0, ylim_top), names.arg = c("Sex","Age","Statins"),
                cex.lab = 1.3, cex.names = 1.3,
                ylab = paste0("Average Portability in ", Ancs[anc]," population"))
  arrows(x0 = bp, y0 = lower, x1 = bp, y1 = upper, angle = 90, code = 3, length = 0.05)
  for (ctx_i in seq_along(context_list)) {
    ctx    <- context_list[ctx_i]
    ci_top <- max(upper[, ctx_i], na.rm=TRUE)
    for (comp in port_comps) {
      x_l   <- bp[comp$i, ctx_i]
      x_r   <- bp[comp$j, ctx_i]
      brk_y <- ci_top + comp$lvl * brk_step - 0.03
      segments(x_l, brk_y, x_r, brk_y)
      segments(x_l, brk_y, x_l, brk_y - tick_len)
      segments(x_r, brk_y, x_r, brk_y - tick_len)
      text(mean(c(x_l, x_r)), brk_y + text_off + 0.005, fmt_p(pval_port[ctx, anc, comp$pv]), cex = 1)
    }
  }
  if (anc == valid_pops[2]){
    legend("topright", legend = c("pgs","ampPGS","PGSC"), bty = "n",
           fill = base_cols3[c('pgs','ampPGS','PGSC_v_pgs')])
    legend("topleft", legend = "(A)", bty="n", inset=c(-0.08,-0.08), xpd=TRUE)
  }else{
    legend("topleft", legend = "(B)", bty="n", inset=c(-0.08,-0.08), xpd=TRUE)
  }
}
dev.off()

# avg_perR2diff 
png("avg_perR2diff.png", width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1,2), mar = c(4,5,2,1))
for (anc in valid_pops[2:3]) {
  means    <- t(avg_R2diff_ivw[, anc, , 'delta'])
  upper    <- t(avg_R2diff_ivw[, anc, , 'delta'] + 1.96 * avg_R2diff_ivw[, anc, , 'se'])
  lower    <- t(avg_R2diff_ivw[, anc, , 'delta'] - 1.96 * avg_R2diff_ivw[, anc, , 'se'])
  ylim_use <- range(lower, upper) + c(-0.5, 1)
  bp <- barplot(means, beside = TRUE, col = c(base_cols3['ampPGS'], base_cols3['PGSC_v_pgs']),
                las = 1, ylim = ylim_use, names.arg = c("Sex","Age","Statins"),
                cex.lab = 1.3, cex.names = 1.3,
                ylab = paste0("Average R² %Change in ", Ancs[anc], " population"))
  arrows(x0 = bp, y0 = lower, x1 = bp, y1 = upper, angle = 90, code = 3, length = 0.05)
  tick_len <- 0.15
  for (ctx_i in seq_along(context_list)) {
    x_amp  <- bp[1, ctx_i]
    x_pgsc <- bp[2, ctx_i]
    brk_y  <- max(upper[, ctx_i], na.rm=TRUE) + 0.3
    segments(x_amp,  brk_y, x_pgsc, brk_y)
    segments(x_amp,  brk_y, x_amp,  brk_y - tick_len)
    segments(x_pgsc, brk_y, x_pgsc, brk_y - tick_len)
    text(mean(c(x_amp, x_pgsc)), brk_y + 0.3, fmt_p(pval_R2diff_ivw[context_list[ctx_i], anc]), cex = 1)
  }
  if (anc == valid_pops[2]){
    legend("topleft", legend = "(A)", bty="n", inset=c(-0.08,-0.1), xpd=TRUE)
  }else{
    legend("right", legend = c("ampPGS","PGSC"), bty = "n",
           fill = c(base_cols3['ampPGS'], base_cols3['PGSC_v_pgs']))
    legend("topleft", legend = "(B)", bty="n", inset=c(-0.08,-0.1), xpd=TRUE)
  }
}
dev.off()
