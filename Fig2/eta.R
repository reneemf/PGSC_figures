rm( list=ls() )
library(dplyr)

source("defaults.R")
source("functions.R")
source("build_sim.R")

### simulate changes to heteroskedasticity
args <- commandArgs(trailingOnly = TRUE)
j <- as.integer(args[1])

gxe_values <- c(0, 0.05)
eta_values <- 10^seq(-2, 2, length = tuner+1)

# per-iteration output file
dir.create('Rdata/eta/', showWarnings=FALSE)
for(n in seq_along(gxe_values)){
  gxe_fixed <- gxe_values[n]
  iter_savefile <- paste0('Rdata/eta/eta_gxe_', gxe_fixed, '_iter_', j, '.Rdata')
  if(!file.exists(iter_savefile)){
    r2_iter <- array(NA, dim=c(length(eta_values),length(pgs_methods)),
                  dimnames=list(eta_values,pgs_methods))
    r2_diff_iter <- array(NA, dim=c(length(eta_values),length(pgs_methods[-1])),
                          dimnames=list(eta_values,pgs_methods[-1]))

    for(i in seq_along(eta_values)){
      temp <- run_sim(N, S, S_caus, h2_add, gxe_fixed,
                                 h2_coord, eta_values[i], alpha, lambda_boxcox)
      r2_iter[i,] <- temp[1:3]
      r2_diff_iter[i,] <- ((temp[2:3] - temp[1])/temp[1])*100
    }

    save( r2_iter, r2_diff_iter, j, eta_values, pgs_methods, file=iter_savefile)
  }
}
