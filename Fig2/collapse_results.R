rm(list = ls())
## Collapse per-iteration Rdata files into single arrays per sweep
## Required modules: module load gcc/12.1.0   module load R/4.3.1

setwd(file.path(Sys.getenv("PGSC_HOME", unset = "."), "Fig2"))
source("defaults.R")

# Parameter grids (must match sweep scripts)
gxe_values   <- c(0, 0.05)
alpha_values <- c(0.05, 0.1, 0.2, 0.5)
LB_values    <- c(0, 1, 2)

# Standard collapse: stacks r2_iter (n_grid x n_pgs) across iterations
collapse_family <- function(dir, pattern, grid_var, out_basename) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) != expected_n_iter) {
    message("Skipping ", dir, ": missing iterations (found ", length(files),
            ", expected ", expected_n_iter, ")")
    return(invisible(NULL))
  }
  files <- sort(files)

  load(files[1])  # load first file to get dimensions
  grid_vals <- get(grid_var)
  n_grid    <- length(grid_vals)
  n_pgs     <- length(pgs_methods)
  n_iter    <- length(files)

  r2_all <- array(NA, dim = c(n_grid, n_pgs, n_iter),
                  dimnames = list(grid_vals, pgs_methods, iter = seq_len(n_iter)))
  r2_diff_all <- array(NA, dim = c(n_grid, n_pgs - 1L, n_iter),
                       dimnames = list(grid_vals, pgs_methods[-1L], iter = seq_len(n_iter)))
  iter_ids <- integer(n_iter)

  for (k in seq_along(files)) {
    load(files[k])
    r2_all[,, k]      <- r2_iter
    r2_diff_all[,, k] <- r2_diff_iter
    iter_ids[k]       <- if (exists("j")) j else k
  }

  dimnames(r2_all)[[3]]      <- iter_ids
  dimnames(r2_diff_all)[[3]] <- iter_ids

  out_file <- file.path(dir, paste0(out_basename, "_collapsed.Rdata"))
  save(r2_all, r2_diff_all, grid_vals, pgs_methods, iter_ids, file = out_file)
  message("Wrote: ", out_file)

  removed <- file.remove(files)
  message("Removed ", sum(removed), " per-iteration file(s)")
}

# c_prob collapse: also stacks r2_strat_iter (stratified C-group R2s)
collapse_c_prob <- function(dir, pattern, grid_var, out_basename) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) != expected_n_iter) {
    message("Skipping ", dir, ": missing iterations (found ", length(files),
            ", expected ", expected_n_iter, ")")
    return(invisible(NULL))
  }
  files <- sort(files)

  load(files[1])
  grid_vals     <- get(grid_var)
  strat_methods <- c("pgs0", "pgs1", "PGSC0", "PGSC1")
  n_grid        <- length(grid_vals)
  n_iter        <- length(files)

  r2_all       <- array(NA, dim = c(n_grid, length(pgs_methods), n_iter),
                        dimnames = list(grid_vals, pgs_methods, seq_len(n_iter)))
  r2_strat_all <- array(NA, dim = c(n_grid, 4, n_iter),
                        dimnames = list(grid_vals, strat_methods, seq_len(n_iter)))
  r2_diff_all  <- array(NA, dim = c(n_grid, length(pgs_methods) - 1L, n_iter),
                        dimnames = list(grid_vals, pgs_methods[-1L], seq_len(n_iter)))
  iter_ids <- integer(n_iter)

  for (k in seq_along(files)) {
    load(files[k])
    r2_all[,, k]       <- r2_iter
    r2_strat_all[,, k] <- r2_strat_iter
    r2_diff_all[,, k]  <- r2_diff_iter
    iter_ids[k]        <- if (exists("j")) j else k
  }

  dimnames(r2_all)[[3]]       <- iter_ids
  dimnames(r2_strat_all)[[3]] <- iter_ids
  dimnames(r2_diff_all)[[3]]  <- iter_ids

  out_file <- file.path(dir, paste0(out_basename, "_collapsed.Rdata"))
  save(r2_all, r2_strat_all, r2_diff_all, grid_vals, pgs_methods, iter_ids,
       file = out_file)
  message("Wrote: ", out_file)

  removed <- file.remove(files)
  message("Removed ", sum(removed), " per-iteration file(s)")
}

# Single-file sweeps (no GxC-background split)
collapse_family("Rdata/h2_gxe",   pattern = "^h2_gxe_iter_.*\\.Rdata$",
                grid_var = "h2_gxe_values",   out_basename = "h2_gxe")
collapse_family("Rdata/h2_coord", pattern = "^h2_coord_iter_.*\\.Rdata$",
                grid_var = "h2_coord_values", out_basename = "h2_coord")

# Sweeps across GxC backgrounds
for (gxe_fixed in gxe_values) {

  collapse_family("Rdata/eta",
                  pattern  = paste0("^eta_gxe_",  gxe_fixed, "_iter_.*\\.Rdata$"),
                  grid_var = "eta_values",
                  out_basename = paste0("eta_", gxe_fixed))

  collapse_family("Rdata/h2",
                  pattern  = paste0("^h2_gxe_",   gxe_fixed, "_iter_.*\\.Rdata$"),
                  grid_var = "h2_values",
                  out_basename = paste0("h2_", gxe_fixed))

  collapse_family("Rdata/poly",
                  pattern  = paste0("^poly_gxe_", gxe_fixed, "_iter_.*\\.Rdata$"),
                  grid_var = "poly_values",
                  out_basename = paste0("poly_", gxe_fixed))

  collapse_family("Rdata/pop",
                  pattern  = paste0("^pop_gxe_",  gxe_fixed, "_iter_.*\\.Rdata$"),
                  grid_var = "N_values",
                  out_basename = paste0("pop_", gxe_fixed))

  for (alpha_fixed in alpha_values) {
    collapse_family("Rdata/lam_box",
                    pattern  = paste0("^lam_box_gxe_", gxe_fixed, "_alpha_",
                                      alpha_fixed, "_iter_.*\\.Rdata$"),
                    grid_var = "LB_values",
                    out_basename = paste0("lam_box_", gxe_fixed, "_alpha_", alpha_fixed))
  }

  for (LB_fixed in LB_values) {
    collapse_family("Rdata/alpha",
                    pattern  = paste0("^alpha_gxe_", gxe_fixed, "_LB_",
                                      LB_fixed, "_iter_.*\\.Rdata$"),
                    grid_var = "alpha_values",
                    out_basename = paste0("alpha_", gxe_fixed, "_LB_", LB_fixed))
  }

  collapse_c_prob("Rdata/c_prob",
                  pattern  = paste0("^c_prob_gxe_", gxe_fixed, "_iter_.*\\.Rdata$"),
                  grid_var = "c_values",
                  out_basename = paste0("c_prob_", gxe_fixed))
}
