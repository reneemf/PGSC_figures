rm( list=ls() )
library(dplyr)

source("defaults.R")
source("functions.R")
source("build_sim.R")

### simulate effects of C effect across the log & sq scales

args <- commandArgs(trailingOnly = TRUE)
j <- as.integer(args[1])

gxe_values <- c(0, 0.05)
alpha_values <- seq(0,5,length = tuner)
LB_values <- c(0,1,2)

# per-iteration output file
dir.create('Rdata/alpha/', showWarnings=FALSE)
for(m in seq_along(gxe_values)){
  gxe_fixed <- gxe_values[m]
  for(n in seq_along(LB_values)){
    LB_fixed <- LB_values[n]
    iter_savefile <- paste0('Rdata/alpha/alpha_gxe_', gxe_fixed,'_LB_',LB_fixed,'_iter_', j, '.Rdata')
    if(!file.exists(iter_savefile)){
      r2_iter <- array(NA, dim=c(length(alpha_values),length(pgs_methods)),
                    dimnames=list(alpha_values,pgs_methods))
      r2_diff_iter <- array(NA, dim=c(length(alpha_values),length(pgs_methods[-1])),
                            dimnames=list(alpha_values,pgs_methods[-1]))

      for(i in seq_along(alpha_values)){
        temp <- run_sim(N, S, S_caus, h2_add, gxe_fixed,
                                   h2_coord, eta, alpha_values[i], LB_fixed)
        r2_iter[i,] <- temp[1:3]
        r2_diff_iter[i,] <- ((temp[2:3] - temp[1])/temp[1])*100
      }

      save( r2_iter, r2_diff_iter, j, alpha_values, pgs_methods, file=iter_savefile)
    }
  }
}
