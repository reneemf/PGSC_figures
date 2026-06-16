rm(list = ls())
setwd(file.path(Sys.getenv("PGSC_HOME", unset = "."), "Fig2"))
source("defaults.R")

load_collapsed <- function(path) {
  if (!file.exists(path)) return(NULL)
  e <- new.env()
  load(path, envir = e)
  as.list(e)
}

# Compute mean R2 and mean %Change from a collapsed r2_all array.
compute_summary <- function(d) {
  if (is.null(d)) return(NULL)
  mean_r2 <- apply(d$r2_all, 1:2, mean, na.rm = TRUE)
  mean_diff <- apply(d$r2_all[, pgs_methods[2:3], , drop = FALSE] - d$r2_all[, c(1,1), , drop = FALSE],
                     1:2, mean, na.rm = TRUE)
  mean_per <- (mean_diff / mean_r2[, "pgs"]) * 100
  # mean_per <- apply(d$r2_diff_all, 1:2, mean, na.rm = TRUE)
  list(x_range = d$grid_vals, mean_r2 = mean_r2, mean_per = mean_per, n_iter = dim(d$r2_all)[3])
}

# Plot Mean R2 %Change for ampPGS and PGSC.
plot_per_diff <- function(dat, xlab, ylab, methods = pgs_methods[2:3], type = "R2_diff"){
  if(type == "R2") y_vals <- dat$mean_r2 else y_vals <- dat$mean_per
  plot(dat$x_range, y_vals[, methods[1]], type = "n",
       xlab = xlab, ylab = ylab,
       xlim = range(dat$x_range), ylim = range(y_vals[, methods]),
       cex.main = 2.2, cex.lab  = 2.5, cex.axis = 1.8)
  abline(h = 0, col = "black", lty = 2)

  for (m in methods){
    lines(dat$x_range, y_vals[, m], lty = 1, lwd = 3, col = base_cols[m])
  }
}

add_std_legend <- function(pos = "topleft", methods = pgs_methods[2:3], cex = 1.8){
  legend(pos, legend = methods, col = base_cols[methods],
         lty = 1, lwd = 3, bty = "n", cex = cex)
}

add_param_note <- function(...){
  legend("right", bty = "n", cex = 1.8,
         legend = c(list(...))) # bquote(S[caus] ~ ":" ~ .(S_caus))
}

dir.create("figs", showWarnings = FALSE)

# Two sets by GxC background
# gxe0:  h2_coord = h2_uncoord = 0
# gxenz: h2_coord = h2_uncoord = 0.05*h2_add

sets <- list(
  list(label   = "gxe0",
       h2      = "Rdata/h2/h2_gxe0_collapsed.Rdata",
       eta     = "Rdata/eta/eta_gxe0_collapsed.Rdata",
       poly    = "Rdata/poly/poly_gxe0_collapsed.Rdata",
       pop     = "Rdata/pop/pop_gxe0_collapsed.Rdata",
       lam_box = "Rdata/lam_box/lam_box_gxe0_collapsed.Rdata"),
  list(label   = "gxenz",
       h2      = "Rdata/h2/h2_gxenz_collapsed.Rdata",
       eta     = "Rdata/eta/eta_gxenz_collapsed.Rdata",
       poly    = "Rdata/poly/poly_gxenz_collapsed.Rdata",
       pop     = "Rdata/pop/pop_gxenz_collapsed.Rdata",
       lam_box = "Rdata/lam_box/lam_box_gxenz_collapsed.Rdata")
)

# h2_uncoord and h2_coord are run once w 0 gxe
h2uncoord_file     <- "Rdata/h2_gxe/h2_gxe_collapsed.Rdata"
h2coord_file   <- "Rdata/h2_coord/h2_coord_collapsed.Rdata"
d_h2uncoord     <- compute_summary(load_collapsed(h2uncoord_file))
d_h2coord   <- compute_summary(load_collapsed(h2coord_file))

# h2_coord ratio x-axis: proportion of total GxC that is coordinated.
h2coord_ratio_x <- d_h2coord$x_range / max(d_h2coord$x_range)
n_lab <- d_h2uncoord$n_iter

# Panels A & B builder (h2xgxe, prop coord, eta w 0 gxe)
draw_panel_AB <- function(out_path, ylab, d_h2uncoord_in, d_h2coord_in, methods,
                         ratio_x_in, d_eta_in, type = "R2_diff") {
  png(out_path, width = 24, height = 8, units = "in", res = 300)
  par(mfrow = c(1, 3), mar = c(5, 8, 4, 1), mgp = c(4, 0.8, 0),
      cex.lab = 1.5, cex.axis = 1.3)
  # 1: h2_uncoord
    plot_per_diff(d_h2uncoord_in, methods,
                  xlab = expression("Locus-specific GxC heritability" ~ (h[uncoord]^2)),
                  ylab = ylab, type = type)
    add_std_legend("topleft", methods)
    add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[coord]^2 ~ ":" ~ .(h2_coord)))

  # 2: h2_coord/(h2_coord+h2_uncoord) — total GxC fixed at 0.1*h2_add
    d_ratio           <- d_h2coord_in
    d_ratio$x_range <- ratio_x_in
    plot_per_diff(d_ratio, methods,
                  xlab = expression("Proportion coordinated GxC heritability" ~ (h[coord]^2 / h[GxC]^2)),
                  ylab = " ", type = type)
    add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[GxC]^2 ~ "total:" ~ .(h2_gxe)))

  # 3: eta / heteroskedasticity
  plot_per_diff(d_eta_in, methods,
                xlab = expression("Heteroskedasticity" ~ (eta)),
                ylab = " ", type = type)
  add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[GxC]^2 ~ "total:" ~ .(h2_uncoord)))

  dev.off()
}

for (s in sets) {
  d_h2        <- compute_summary(load_collapsed(s$h2))
  d_eta       <- compute_summary(load_collapsed(s$eta))
  d_poly      <- compute_summary(load_collapsed(s$poly))
  d_N         <- compute_summary(load_collapsed(s$pop))
  d_lam_box   <- compute_summary(load_collapsed(s$lam_box))

  if(s$label == "gxe0"){
    total_gxe <- h2_coord

    # Panel A: 1x3 row (h2_uncoord, h2_coord, eta)
    draw_panel_AB(paste0("figs/panel_A_", s$label, ".png"),
                 methods = pgs_methods[2:3],
                 ylab = paste0("Mean R\u00b2 %Change (over ", n_lab, " iterations)"),
                 d_h2uncoord, d_h2coord, h2coord_ratio_x, d_eta)
    # Panel B: 1x3 row (h2_uncoord, h2_coord, eta) R2s
    draw_panel_AB(paste0("figs/panel_B_", s$label, ".png"),
                 methods = pgs_methods[1:3],
                 ylab = paste0("Mean R\u00b2 (over ", n_lab, " iterations)"),
                 d_h2uncoord, d_h2coord, h2coord_ratio_x, d_eta, type = "R2")
  }else{
    total_gxe <- h2_gxe
  }
  # Panel C: 1x4 row (h2, lam_box, poly, pop)
  png(paste0("figs/panel_C_", s$label, ".png"),
      width = 32, height = 8, units = "in", res = 300)
  par(mfrow = c(1, 4), mar = c(5, 8, 4, 1), mgp = c(4, 0.8, 0),
      cex.lab = 1.5, cex.axis = 1.3)

  # 1: h2
  plot_per_diff(d_h2, xlab = expression("Additive heritability" ~ (h^2)),
                ylab = paste0("Mean R\u00b2 %Change (over ", n_lab, " iterations)"))
  add_std_legend("topright")
  add_param_note(bquote(h[GxC]^2 ~ "total:" ~ .(total_gxe)))

  # 2: lambda_boxcox
  plot_per_diff(d_lam_box, xlab = expression("Scaling effect" ~ (lambda[BoxCox])), ylab = " ")
  add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[GxC]^2 ~ "total:" ~ .(total_gxe)))

  # 3: S_caus
  plot_per_diff(d_poly, xlab = expression("Polygenicity" ~ (S[caus])), ylab = " ")
  add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[GxC]^2 ~ "total:" ~ .(total_gxe)))

  # 4: N
  plot_per_diff(d_N, xlab = "Sample size (N)", ylab = " ")
  add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[GxC]^2 ~ "total:" ~ .(total_gxe)))

  dev.off()

  # Panel D: 1x4 row (h2, lam_box, poly, pop) R2s
  png(paste0("figs/panel_D_", s$label, ".png"),
      width = 32, height = 8, units = "in", res = 300)
  par(mfrow = c(1, 4), mar = c(5, 8, 4, 1), mgp = c(4, 0.8, 0),
      cex.lab = 1.5, cex.axis = 1.3)

  # 1: h2
  plot_per_diff(d_h2, xlab = expression("Additive heritability" ~ (h^2)),
                ylab = paste0("Mean R\u00b2 (over ", n_lab, " iterations)"),
                methods = pgs_methods[1:3], type = "R2")
  add_std_legend("topleft", methods = pgs_methods[1:3])
  add_param_note(bquote(h[GxC]^2 ~ "total:" ~ .(total_gxe)))

  # 2: lambda_boxcox
  plot_per_diff(d_lam_box, xlab = expression("Scaling effect" ~ (lambda[BoxCox])),
                ylab = " ", methods = pgs_methods[1:3], type = "R2")
  add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[GxC]^2 ~ "total:" ~ .(total_gxe)))

  # 3: S_caus
  plot_per_diff(d_poly, xlab = expression("Polygenicity" ~ (S[caus])),
                ylab = " ", methods = pgs_methods[1:3], type = "R2")
  add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[GxC]^2 ~ "total:" ~ .(total_gxe)))

  # 4: N
  plot_per_diff(d_N, xlab = "Sample size (N)",
                ylab = " ", methods = pgs_methods[1:3], type = "R2")
  add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[GxC]^2 ~ "total:" ~ .(total_gxe)))

  dev.off()
}


