rm( list=ls() )
library(dplyr)

source("defaults.R")
source("functions.R")
source("build_sim.R")

### simulate changes to proportion C
args <- commandArgs(trailingOnly = TRUE)
j <- as.integer(args[1])

gxe_values <- c(0, 0.05)
c_values <- seq(0.1, 0.9, length = tuner)

# per-iteration output file
dir.create('Rdata/c_prob/', showWarnings=FALSE)
for(n in seq_along(gxe_values)){
  gxe_fixed <- gxe_values[n]
  iter_savefile <- paste0('Rdata/c_prob/c_prob_gxe_', gxe_fixed, '_iter_', j, '.Rdata')
  if(!file.exists(iter_savefile)){
    r2_iter <- array(NA, dim=c(length(c_values), length(pgs_methods)),
                     dimnames=list(c_values, pgs_methods))
    r2_strat_iter <- array(NA, dim=c(length(c_values), 4),
                           dimnames=list(c_values, c("pgs0","pgs1","PGSC0","PGSC1")))
    r2_diff_iter <- array(NA, dim=c(length(c_values),length(pgs_methods[-1])),
                          dimnames=list(c_values,pgs_methods[-1]))

    for(i in seq_along(c_values)){
      temp <- run_sim(N, S, S_caus, h2_add, gxe_fixed,
                                 h2_coord, eta, alpha, lambda_boxcox, c_values[i])
      r2_iter[i,]       <-   temp[1:3]
      r2_strat_iter[i,] <-   temp[4:7]
      r2_diff_iter[i,]  <- ((temp[2:3] - temp[1])/temp[1])*100
    }
    save(r2_iter, r2_strat_iter, r2_diff_iter, j, c_values, pgs_methods, file=iter_savefile)
  }
}
