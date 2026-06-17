rm(list = ls())
setwd("/Users/reneefonseca/Documents/UChicago/Dahl/Figures/Fig2/")
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
  mean_diff <- apply(d$r2_all[, pgs_methods[2:5], , drop = FALSE] - d$r2_all[, c(1,1,1,1), , drop = FALSE],
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

# empty plot when data is missing
plot_or_empty <- function(dat, xlab, ylab, methods = pgs_methods[2:3], type = "R2_diff"){
  if (is.null(dat)) {
    plot(0, 0, type = "n", xlab = xlab, ylab = ylab,
         xaxt = "n", yaxt = "n", cex.lab = 2.5)
    text(0, 0, " ", cex = 2.5, col = "grey60")
    return(invisible())
  }
  plot_per_diff(dat, xlab = xlab, ylab = ylab, methods = methods, type = type)
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

# Two sets by GxE background
# gxe0:  h2_coord = h2_uncoord = 0
# gxenz: h2_coord = h2_uncoord = 0.05

sets <- list(
  list(label   = "gxe0",
       alpha0  = "Rdata/alpha/alpha_0_LB_0_collapsed.Rdata",
       alpha1  = "Rdata/alpha/alpha_0_LB_1_collapsed.Rdata",
       alpha2  = "Rdata/alpha/alpha_0_LB_2_collapsed.Rdata",
       c_prob  = "Rdata/c_prob/c_prob_0_collapsed.Rdata",
       h2      = "Rdata/h2/h2_0_collapsed.Rdata",
       eta     = "Rdata/eta/eta_0_collapsed.Rdata",
       poly    = "Rdata/poly/poly_0_collapsed.Rdata",
       pop     = "Rdata/pop/pop_0_collapsed.Rdata",
       lam_box = "Rdata/lam_box/lam_box_0_collapsed.Rdata"),
  list(label   = "gxenz",
       alpha0   = "Rdata/alpha/alpha_0.05_LB_0_collapsed.Rdata",
       alpha1   = "Rdata/alpha/alpha_0.05_LB_1_collapsed.Rdata",
       alpha2   = "Rdata/alpha/alpha_0.05_LB_2_collapsed.Rdata",
       c_prob  = "Rdata/c_prob/c_prob_0.05_collapsed.Rdata",
       h2      = "Rdata/h2/h2_0.05_collapsed.Rdata",
       eta     = "Rdata/eta/eta_0.05_collapsed.Rdata",
       poly    = "Rdata/poly/poly_0.05_collapsed.Rdata",
       pop     = "Rdata/pop/pop_0.05_collapsed.Rdata",
       lam_box = "Rdata/lam_box/lam_box_0.05_collapsed.Rdata")
)

# h2_uncoord and h2_coord are run once w 0 gxe
h2uncoord_file  <- "Rdata/h2_gxe/h2_gxe_collapsed.Rdata"
h2coord_file    <- "Rdata/h2_coord/h2_coord_collapsed.Rdata"
d_h2uncoord     <- compute_summary(load_collapsed(h2uncoord_file))
d_h2coord       <- compute_summary(load_collapsed(h2coord_file))

# h2_coord ratio x-axis: proportion of total GxE that is coordinated.
h2coord_ratio_x <- if (!is.null(d_h2coord)) d_h2coord$x_range / max(d_h2coord$x_range) else NULL
n_lab <- if (!is.null(d_h2uncoord)) d_h2uncoord$n_iter else expected_n_iter

# Panels A & B builder (h2_uncoord, prop coord, eta -- all w/ 0 gxe)
draw_panel_AB <- function(out_path, ylab, d_h2uncoord_in, d_h2coord_in, methods,
                         ratio_x_in, d_eta_in, type = "R2_diff") {
  png(out_path, width = 24, height = 8, units = "in", res = 300)
  par(mfrow = c(1, 3), mar = c(5, 8, 4, 1), mgp = c(4, 0.8, 0),
      cex.lab = 1.5, cex.axis = 1.3)
  # 1: h2_uncoord
    plot_or_empty(d_h2uncoord_in, methods,
                  xlab = expression("Locus-specific GxE heritability" ~ (h[uncoord]^2)),
                  ylab = ylab, type = type)
    add_std_legend("topleft", methods)
    add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[coord]^2 ~ ":" ~ .(h2_coord)))

  # 2: h2_coord/(h2_coord+h2_uncoord) 
    d_ratio         <- d_h2coord_in
    if (!is.null(d_ratio)) d_ratio$x_range <- ratio_x_in
    plot_or_empty(d_ratio, methods,
                  xlab = expression("Proportion coordinated GxE heritability" ~ (h[coord]^2 / h[GxE]^2)),
                  ylab = " ", type = type)
    add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[GxE]^2 ~ ":" ~ .(h2_gxe)))

  # 3: eta / heteroskedasticity
  plot_or_empty(d_eta_in, methods,
                xlab = expression("Heteroskedasticity" ~ (eta)),
                ylab = " ", type = type)
  add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[GxE]^2 ~ ":" ~ .(h2_uncoord)))

  dev.off()
}

# Panels C/D/E/F builder: 2x4 grid sweeping
#   h2, polygenicity, lambda_boxcox, c_prob, alpha LB_0, alpha LB_1, alpha LB_2.
# `d` is a named list of compute_summary() outputs (any may be NULL/pending).
draw_panel_grid <- function(out_path, d, ylab, methods, total_gxe, type = "R2_diff") {
  items <- list(
    list(dat = d$h2,      xlab = expression("Additive heritability" ~ (h^2))),
    list(dat = d$poly,    xlab = expression("Polygenicity" ~ (S[caus]))),
    list(dat = d$lam_box, xlab = expression("Scaling effect" ~ (lambda[BoxCox]))),
    list(dat = d$c_prob,  xlab = expression("Distribution of C" ~ (c[prob]))),
    list(dat = d$alpha0,  xlab = expression("Main effect of C" ~ (alpha)), note = "Log scale"),
    list(dat = d$alpha1,  xlab = expression("Main effect of C" ~ (alpha)), note = "True scale"),
    list(dat = d$alpha2,  xlab = expression("Main effect of C" ~ (alpha)), note = "Squared scale")
  )

  png(out_path, width = 32, height = 16, units = "in", res = 300)
  par(mfrow = c(2, 4), mar = c(5, 8, 4, 1), mgp = c(4, 0.8, 0),
      cex.lab = 1.5, cex.axis = 1.3)
  for (i in seq_along(items)) {
    it <- items[[i]]
    yl <- if (i %in% c(1, 5)) ylab else " "
    plot_or_empty(it$dat, xlab = it$xlab, ylab = yl, methods = methods, type = type)
    if (i == 1) add_std_legend("topright", methods)
    notes <- list(bquote(h[GxE]^2 ~ "total:" ~ .(total_gxe)))
    if (!is.null(it$note)) notes <- c(notes, list(it$note))
    do.call(add_param_note, notes)
  }
  dev.off()
}

for (s in sets) {
  d <- list(
    alpha0  = compute_summary(load_collapsed(s$alpha0)),
    alpha1  = compute_summary(load_collapsed(s$alpha1)),
    alpha2  = compute_summary(load_collapsed(s$alpha2)),
    c_prob  = compute_summary(load_collapsed(s$c_prob)),
    h2      = compute_summary(load_collapsed(s$h2)),
    eta     = compute_summary(load_collapsed(s$eta)),
    poly    = compute_summary(load_collapsed(s$poly)),
    lam_box = compute_summary(load_collapsed(s$lam_box))
  )

  if (s$label == "gxe0") {
    total_gxe <- 0

    ##### Panel A: 1x3 row (h2_uncoord, prop coord, eta) -- R2 %Change
    draw_panel_AB(paste0("figs/panel_A_", s$label, ".png"),
                 methods = pgs_methods[2:3],
                 ylab = paste0("Mean R² %Change (over ", n_lab, " iterations)"),
                 d_h2uncoord, d_h2coord, h2coord_ratio_x, d$eta)
    ##### Panel B: 1x3 row (h2_uncoord, prop coord, eta) -- R2
    draw_panel_AB(paste0("figs/panel_B_", s$label, ".png"),
                 methods = pgs_methods[1:3],
                 ylab = paste0("Mean R² (over ", n_lab, " iterations)"),
                 d_h2uncoord, d_h2coord, h2coord_ratio_x, d$eta, type = "R2")

    ##### Panel C: 2x4 grid -- R2 %Change (ampPGS, PGSC)
    draw_panel_grid(paste0("figs/panel_C_", s$label, ".png"), d,
                    ylab = paste0("Mean R² %Change (over ", n_lab, " iterations)"),
                    methods = pgs_methods[2:3], total_gxe = total_gxe)
    ##### Panel D: 2x4 grid -- R2 (pgs, ampPGS, PGSC)
    draw_panel_grid(paste0("figs/panel_D_", s$label, ".png"), d,
                    ylab = paste0("Mean R² (over ", n_lab, " iterations)"),
                    methods = pgs_methods[1:3], total_gxe = total_gxe, type = "R2")
  } else {
    total_gxe <- 0.05

    ##### Panel E: 2x4 grid -- R2 %Change (ampPGS, PGSC)
    draw_panel_grid(paste0("figs/panel_E_", s$label, ".png"), d,
                    ylab = paste0("Mean R² %Change (over ", n_lab, " iterations)"),
                    methods = pgs_methods[2:3], total_gxe = total_gxe)
    ##### Panel F: 2x4 grid -- R2 (pgs, ampPGS, PGSC)
    draw_panel_grid(paste0("figs/panel_F_", s$label, ".png"), d,
                    ylab = paste0("Mean R² (over ", n_lab, " iterations)"),
                    methods = pgs_methods[1:3], total_gxe = total_gxe, type = "R2")
  }
}
