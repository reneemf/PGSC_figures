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
  n_alt <- length(pgs_methods) - 1L
  mean_diff <- apply(d$r2_all[, pgs_methods[-1], , drop = FALSE] - d$r2_all[, rep(1L, n_alt), , drop = FALSE],
                     1:2, mean, na.rm = TRUE)
  mean_per <- (mean_diff / mean_r2[, "pgs"]) * 100
  # mean_per <- apply(d$r2_diff_all, 1:2, mean, na.rm = TRUE)
  list(x_range = d$grid_vals, mean_r2 = mean_r2, mean_per = mean_per, n_iter = dim(d$r2_all)[3])
}

# Plot Mean R2 %Change for ampPGS and PGSC.
# logx = TRUE puts the x-axis on a log scale. logx_offset is added to every x
# value first so a 0 (log(0) = -Inf) still plots -- keeps all points.
plot_per_diff <- function(dat, xlab, ylab, methods = pgs_methods[2:3], type = "R2_diff", ylim = NULL, logx = FALSE, logx_offset = 0){#.001){
  if(type == "R2") y_vals <- dat$mean_r2 else y_vals <- dat$mean_per
  if (is.null(ylim)) ylim <- range(y_vals[, methods])
  x <- if (logx) dat$x_range + logx_offset else dat$x_range
  plot(x, y_vals[, methods[1]], type = "n",
       xlab = xlab, ylab = ylab,
       xlim = range(x), ylim = ylim,
       log = if (logx) "x" else "",
       cex.main = 3, cex.lab  = 3, cex.axis = 2.5)
  abline(h = 0, col = "black", lty = 2)

  for (m in methods){
    lines(x, y_vals[, m], lty = 1, lwd = 4, col = base_cols[m])
  }
}

# Plot wrapper: draws a "data pending" placeholder when a sweep hasn't been
# generated yet (compute_summary returned NULL), so the panel still renders.
plot_or_empty <- function(dat, xlab, ylab, methods = pgs_methods[2:3], type = "R2_diff", ylim = NULL, logx = FALSE, logx_offset = 0.05){
  if (is.null(dat)) {
    plot(0, 0, type = "n", xlab = xlab, ylab = ylab,
         xaxt = "n", yaxt = "n", cex.lab = 3)
    text(0, 0, "data pending", cex = 3, col = "grey60")
    return(invisible())
  }
  plot_per_diff(dat, xlab = xlab, ylab = ylab, methods = methods, type = type, ylim = ylim, logx = logx, logx_offset = logx_offset)
}

add_std_legend <- function(pos = "topleft", methods = pgs_methods[2:3], cex = 3){
  legend(pos, legend = methods, col = base_cols[methods],
         lty = 1, lwd = 4, bty = "n", cex = cex)
}

add_param_note <- function(..., pos = "right", cex = 2, inset = 0){
  legend(pos, bty = "n", cex = cex, inset = inset,
         legend = c(list(...))) # bquote(S[caus] ~ ":" ~ .(S_caus))
}

add_panel_letter <- function(i){
  mtext(paste0("(", letters[i], ")"), side = 3, line = 1, adj = 0, font = 1, cex = 2)
}

dir.create("figs", showWarnings = FALSE)

# Two sets by GxC background
# gxe0:  h2_coord = h2_uncoord = 0
# gxenz: h2_coord = h2_uncoord = 0.05*h2_add

sets <- list(
  list(label   = "gxe0",
       h2      = "Rdata/h2/h2_0_collapsed.Rdata",
       eta     = "Rdata/eta/eta_0_collapsed.Rdata",
       poly    = "Rdata/poly/poly_0_collapsed.Rdata",
       pop     = "Rdata/pop/pop_0_collapsed.Rdata",
       lam_box = "Rdata/lam_box/lam_box_0_alpha_0.1_collapsed.Rdata"),  # scaling effect: sigma2_C = 0.1
  list(label   = "gxenz",
       h2      = "Rdata/h2/h2_0.05_collapsed.Rdata",
       eta     = "Rdata/eta/eta_0.05_collapsed.Rdata",
       poly    = "Rdata/poly/poly_0.05_collapsed.Rdata",
       pop     = "Rdata/pop/pop_0.05_collapsed.Rdata",
       lam_box = "Rdata/lam_box/lam_box_0.05_alpha_0.1_collapsed.Rdata")  # scaling effect: sigma2_C = 0.1
)

# h2_uncoord and h2_coord are run once w 0 gxe
h2uncoord_file     <- "Rdata/h2_gxe/h2_gxe_collapsed.Rdata"
h2coord_file   <- "Rdata/h2_coord/h2_coord_collapsed.Rdata"
d_h2uncoord     <- compute_summary(load_collapsed(h2uncoord_file))
d_h2coord   <- compute_summary(load_collapsed(h2coord_file))

# h2_coord ratio x-axis: proportion of total GxE that is coordinated.
h2coord_ratio_x <- d_h2coord$x_range / max(d_h2coord$x_range)
n_lab <- d_h2uncoord$n_iter

# Panels A & B builder (h2_uncoord, prop coord, eta w 0 gxe)
draw_panel_AB <- function(out_path, ylab, d_h2uncoord_in, d_h2coord_in, methods,
                         ratio_x_in, d_eta_in, type = "R2_diff", log_eta = FALSE) {
  png(out_path, width = 24, height = 8, units = "in", res = 300)
  par(mfrow = c(1, 3), mar = c(9, 8.1, 4, 4.5), mgp = c(6, 2, 0),
      cex.lab = 3, cex.axis = 2.5, lwd = 2, tcl = -1)
  # 1: h2_uncoord
    plot_or_empty(d_h2uncoord_in, methods,
                  xlab = "Locus-specific GxC heritability",
                  ylab = ylab, type = type)
    add_panel_letter(1)
    add_std_legend("topleft", methods)
    #add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[amp]^2 ~ ":" ~ .(h2_coord)))

  # 2: h2_coord/(h2_coord+h2_uncoord) -- total GxE fixed at 0.1*h2_add
    d_ratio           <- d_h2coord_in
    if (!is.null(d_ratio)) d_ratio$x_range <- ratio_x_in
    plot_or_empty(d_ratio, methods,
                  xlab = expression("Proportion amplification GxC heritability" ~ ( "%" * amp)),
                  ylab = " ", type = type)
    add_panel_letter(2)
    #add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[GxC]^2 + h[amp]^2 ~ "total:" ~ .(h2_gxe)))

  # 3: eta / heteroskedasticity (optionally log-scaled x-axis, eta shifted +0.001)
  plot_or_empty(d_eta_in, methods,
                xlab = if (log_eta) expression("log(Heteroskedasticity)" ~ (eta))
                       else expression("Heteroskedasticity" ~ (eta)),
                ylab = " ", type = type, logx = log_eta)
  add_panel_letter(3)
  #add_param_note(bquote(h^2 ~ ":" ~ .(h2_add)), bquote(h[GxC]^2 ~ ":" ~ .(h2_uncoord)), bquote(h[amp]^2 ~ ":" ~ .(h2_uncoord)))

  dev.off()
}

# Panels C & D builder: 2x4 grid where each row is the SAME four sweeps but a
# different GxE background -- top row = 0 GxE, bottom row = 0.05 h2_GxE.
# `cols` describes the four columns: which sweep to pull from the per-background
# list (`key`), its x-axis label (`xlab`), and an optional annotation (`note`).
# Drop grid points from a column's data. `col$drop_x` lists x-values to remove,
# matched (with tolerance) against the transformed x (col$xform applied) so they
# can be specified on the same scale the axis is drawn in.
drop_col_points <- function(dat, col) {
  if (is.null(dat) || is.null(col$drop_x)) return(dat)
  fx   <- if (!is.null(col$xform)) col$xform(dat$x_range) else dat$x_range
  keep <- !vapply(fx, function(v) any(abs(v - col$drop_x) < 1e-3), logical(1))
  dat$x_range  <- dat$x_range[keep]
  dat$mean_r2  <- dat$mean_r2[keep, , drop = FALSE]
  dat$mean_per <- dat$mean_per[keep, , drop = FALSE]
  dat
}

draw_panel_2gxe <- function(out_path, d_top, d_bot, cols, ylab, methods, type = "R2_diff",
                            col_groups = NULL, legend_pos = "topright", note_pos = "right",
                            legend_cex = 3, note_cex = 2, note_inset = 0) {
  png(out_path, width = 32, height = 16, units = "in", res = 300)
  par(mfrow = c(2, 4), mar = c(9, 9, 4, 1), mgp = c(6, 2, 0),
      cex.lab = 2.5, cex.axis = 2.5, lwd = 2, tcl = -1)
  rows <- list(list(d = d_top, gxe = 0), list(d = d_bot, gxe = 0.05))

  # Shared y-axis ranges: col_groups assigns each column to a group id; all
  # panels (both rows) sharing a group id get one common ylim spanning them.
  group_ylim <- NULL
  if (!is.null(col_groups)) {
    group_ylim <- lapply(unique(col_groups), function(g) {
      vals <- unlist(lapply(rows, function(r)
        lapply(which(col_groups == g), function(j) {
          dat <- drop_col_points(r$d[[cols[[j]]$key]], cols[[j]])
          if (is.null(dat)) return(NULL)
          y <- if (type == "R2") dat$mean_r2 else dat$mean_per
          y[, methods]
        })))
      if (length(vals)) range(vals, na.rm = TRUE) else NULL
    })
    names(group_ylim) <- as.character(unique(col_groups))
  }

  panel_i <- 0
  for (r in rows) {
    for (j in seq_along(cols)) {
      panel_i <- panel_i + 1
      col <- cols[[j]]
      yl  <- if (j == 1) ylab else " "   # ylab only on the left column of each row
      ylim <- if (!is.null(col_groups)) group_ylim[[as.character(col_groups[j])]] else NULL
      dat <- drop_col_points(r$d[[col$key]], col)
      # optional x-axis transform (e.g. main effect -> fraction of variance)
      if (!is.null(dat) && !is.null(col$xform)) dat$x_range <- col$xform(dat$x_range)
      plot_or_empty(dat, xlab = col$xlab, ylab = yl,
                    methods = methods, type = type, ylim = ylim)
      add_panel_letter(panel_i)
      if (panel_i == 1) add_std_legend(legend_pos, methods, cex = legend_cex)
      notes <- list(bquote(h[GxC]^2 ~ ":" ~ .(r$gxe))) #bquote(h[GxC]^2 ~ ":" ~ .(h2_uncoord)), bquote(h[amp]^2 ~ ":" ~ .(h2_uncoord))
      if (!is.null(col$note)) notes <- c(notes, list(col$note))
      do.call(add_param_note, c(notes, list(pos = note_pos, cex = note_cex, inset = note_inset)))
    }
  }
  dev.off()
}

# Load the eight sweeps for one GxE background into a named list.
load_grid <- function(s, gxe_tok) {
  list(
    h2      = compute_summary(load_collapsed(s$h2)),
    poly    = compute_summary(load_collapsed(s$poly)),
    pop     = compute_summary(load_collapsed(s$pop)),
    lam_box = compute_summary(load_collapsed(s$lam_box)),
    c_prob  = compute_summary(load_collapsed(paste0("Rdata/c_prob/c_prob_", gxe_tok, "_collapsed.Rdata"))),
    alpha0  = compute_summary(load_collapsed(paste0("Rdata/alpha/alpha_", gxe_tok, "_LB_0_collapsed.Rdata"))),
    alpha1  = compute_summary(load_collapsed(paste0("Rdata/alpha/alpha_", gxe_tok, "_LB_1_collapsed.Rdata"))),
    alpha2  = compute_summary(load_collapsed(paste0("Rdata/alpha/alpha_", gxe_tok, "_LB_2_collapsed.Rdata")))
  )
}

d0  <- load_grid(sets[[1]], "0")      # 0 GxE background (top rows)
d05 <- load_grid(sets[[2]], "0.05")   # 0.05 h2_GxE background (bottom rows)
d_eta_gxe0 <- compute_summary(load_collapsed(sets[[1]]$eta))

ylab_diff <- paste0("Mean R² %Change (over ", n_lab, " iterations)")
ylab_r2   <- paste0("Mean R² (over ", n_lab, " iterations)")

##### Panels A & B: 1x3 (h2_uncoord, proportion coordinated, heteroskedasticity)
# A: R2 %Change (ampPGS, PGSC)
draw_panel_AB("figs/panel_A.png", methods = pgs_methods[2:3], ylab = ylab_diff,
              d_h2uncoord, d_h2coord, h2coord_ratio_x, d_eta_gxe0)
# A (log eta): same as A but the heteroskedasticity (3rd) plot has a log x-axis.
draw_panel_AB("figs/panel_A_logeta.png", methods = pgs_methods[2:3], ylab = ylab_diff,
              d_h2uncoord, d_h2coord, h2coord_ratio_x, d_eta_gxe0, log_eta = TRUE)
# B: R2 (pgs, ampPGS, PGSC)
draw_panel_AB("figs/panel_B.png", methods = pgs_methods[1:3], ylab = ylab_r2,
              d_h2uncoord, d_h2coord, h2coord_ratio_x, d_eta_gxe0, type = "R2")

##### Panel C: simulation characteristics -- h2, polygenicity, sample size,
# context imbalance. Top row = 0 GxC, bottom row = 0.05 h2_GxE.
cols_char <- list(
  list(key = "h2",     xlab = expression("Additive heritability" ~ (h^2))),
  list(key = "poly",   xlab = expression("Polygenicity" ~ (S[caus]))),
  list(key = "pop",    xlab = "Sample size (N)"),
  list(key = "c_prob", xlab = expression("Context imbalance" ~ (c[prob])))
)
draw_panel_2gxe("figs/panel_C_diff.png", d0, d05, cols_char,
                ylab = ylab_diff, methods = pgs_methods[2:3],
                legend_pos = "bottomright", legend_cex = 3, note_cex = 2.5)
draw_panel_2gxe("figs/panel_C_R2.png",   d0, d05, cols_char,
                ylab = ylab_r2,   methods = pgs_methods[1:3], type = "R2",
                legend_pos = "bottomright", legend_cex = 3, note_cex = 2.5)

##### Panel D: scaling effect (sigma2_C = 0.1) + main effect of C on the log /
# true / squared scales. Top row = 0 GxC, bottom row = 0.05 h2_GxE.
# Main effect of C shown as fraction of total phenotypic variance explained by
# sex: frac = x^2 / (x^2 + 1), where x is the raw main-effect value.
frac_var <- function(x) x^2 / (x^2 + 1)
xlab_frac <- expression("Variance explained by sex" ~ (alpha^2 / (alpha^2 + 1)))
cols_Ceff <- list(
  list(key = "lam_box", xlab = expression("Scaling effect" ~ (lambda[BoxCox])),
       note = bquote(alpha^2 ~ ":" ~ 0.1)),
  list(key = "alpha0",  xlab = xlab_frac, xform = frac_var, note = "Log scale"),
  list(key = "alpha1",  xlab = xlab_frac, xform = frac_var, note = "True scale"),
  list(key = "alpha2",  xlab = xlab_frac, xform = frac_var, note = "Squared scale")
)
# Shared y-axis ranges across rows: (a,e) together, (b,d,f,h) together, (c,g)
# together -- col_groups assigns cols 1/2/3/4 to groups 1/2/3/2.
# Drop the two crowded high-x points (0.917, 0.952 on the alpha^2/(alpha^2+1)
# axis) from plots c & g (the alpha1 / True-scale column). Applied to both the
# diff and R2 panel D figures so they stay consistent.
cols_Ceff_diff <- cols_Ceff
cols_Ceff_diff[[3]]$drop_x <- c(0.917, 0.952)
draw_panel_2gxe("figs/panel_D_diff.png", d0, d05, cols_Ceff_diff,
                ylab = ylab_diff, methods = pgs_methods[2:3],
                col_groups = c(1, 2, 3, 2),
                legend_pos = "right", note_pos = "top",
                legend_cex = 3, note_cex = 2.5, note_inset = c(0, 0.02))
draw_panel_2gxe("figs/panel_D_R2.png",   d0, d05, cols_Ceff_diff,
                ylab = ylab_r2,   methods = pgs_methods[1:3], type = "R2",
                col_groups = c(1, 2, 3, 2),
                legend_pos = "right", note_pos = "top",
                legend_cex = 3, note_cex = 2.5, note_inset = c(0, 0.02))
