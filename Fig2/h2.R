rm( list=ls() )
library(dplyr)

source("defaults.R")
source("functions.R")
source("build_sim.R")

### simulate changes to additive heritability when total GxC heritability is set to 0

# --- get iteration index from command line ---
args <- commandArgs(trailingOnly = TRUE)
j <- as.integer(args[1])

gxe_values <- c(0, 0.05)
h2_values <- seq(0.01, 0.8, length = tuner)

# per-iteration output file
dir.create('Rdata/h2/', showWarnings=FALSE)
for(n in seq_along(gxe_values)){
 gxe_fixed <- gxe_values[n]
  iter_savefile <- paste0('Rdata/h2/h2_gxe_', gxe_fixed, '_iter_', j, '.Rdata')
  if(!file.exists(iter_savefile)){
    r2_iter <- array(NA, dim=c(length(h2_values),length(pgs_methods)),
                  dimnames=list(h2_values,pgs_methods))
    r2_diff_iter <- array(NA, dim=c(length(h2_values),length(pgs_methods[-1])),
                          dimnames=list(h2_values,pgs_methods[-1]))

    for(i in seq_along(h2_values)){
      temp <- run_sim(N, S, S_caus, h2_values[i], gxe_fixed,
                                     h2_coord, eta, alpha, lambda_boxcox)
      r2_iter[i,] <- temp[1:3]
      r2_diff_iter[i,] <- ((temp[2:3] - temp[1])/temp[1])*100
    }

    save( r2_iter, r2_diff_iter, j, h2_values, pgs_methods, file=iter_savefile )
  }
}
