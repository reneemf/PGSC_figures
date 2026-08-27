rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
load(paste0(rdata_dir, "summ_stats.Rdata"))

setwd(paste0(fig_dir, "Fig5"))
context       <- context_list[1]            # "sex"
clean_context <- context_mapping[[context]]$clean_context

# Loaded from summ_stats.Rdata:
#   biome_r2diff_means, biome_r2diff_ses, biome_r2diff_pval  → avg_perR2diff plot
#   biome_port_means, biome_port_ses, biome_port_pvals       → avg_portability plot
#   base_cols3, alt_pgs3, fmt_p

#Plot 1: avg_perR2diff
upper    <- biome_r2diff_means + 1.96 * biome_r2diff_ses
lower    <- biome_r2diff_means - 1.96 * biome_r2diff_ses
ylim_use <- range(lower, upper) + c(-0.5, 0.75)

png("avg_perR2diff_biome.png", width = 5, height = 5, units = "in", res = 300)
par(mar = c(4, 5, 2, 1))
bp <- barplot(biome_r2diff_means, col = base_cols3[names(biome_r2diff_means)],
              ylim = ylim_use, names.arg = c("ampPGS", "PGSC"),
              ylab = "Average R\u00b2 %Change in MSM population (GxSex)")
arrows(x0 = bp, y0 = lower, x1 = bp, y1 = upper, angle = 90, code = 3, length = 0.05)
brk_y    <- max(upper) + 0.3
tick_len <- 0.15
segments(bp[1], brk_y, bp[2], brk_y)
segments(bp[1], brk_y, bp[1], brk_y - tick_len)
segments(bp[2], brk_y, bp[2], brk_y - tick_len)
text(mean(bp), brk_y + 0.3, fmt_p(biome_r2diff_pval), cex = 0.7, xpd=T)
dev.off()

