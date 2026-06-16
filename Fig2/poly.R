rm( list=ls() )
library(dplyr)

source("defaults.R")
source("functions.R")
source("build_sim.R")

### simulate changes to polygenicity by varying S_caus

args <- commandArgs(trailingOnly = TRUE)
j <- as.integer(args[1])

gxe_values <- c(0, 0.05)
poly_values <- seq(S*0.025, S*0.5, length = tuner)

# per-iteration output file
dir.create('Rdata/poly/', showWarnings=FALSE)
for(n in seq_along(gxe_values)){
  gxe_fixed <- gxe_values[n]
  iter_savefile <- paste0('Rdata/poly/poly_gxe_', gxe_fixed, '_iter_', j, '.Rdata')
  if(!file.exists(iter_savefile)){
    r2_iter <- array(NA, dim=c(length(poly_values),length(pgs_methods)),
                  dimnames=list(poly_values,pgs_methods))
    r2_diff_iter <- array(NA, dim=c(length(poly_values),length(pgs_methods[-1])),
                          dimnames=list(poly_values,pgs_methods[-1]))

    for(i in seq_along(poly_values)){
      temp <- unlist(run_sim(N, S, poly_values[i], h2_add, gxe_fixed,
                                 h2_coord, eta, alpha, lambda_boxcox))
      r2_iter[i,] <- temp[1:3]
      r2_diff_iter[i,] <- ((temp[2:3] - temp[1])/temp[1])*100
    }

    save( r2_iter, r2_diff_iter, j, poly_values, pgs_methods, file=iter_savefile)
  }
}
