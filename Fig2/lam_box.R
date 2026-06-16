rm( list=ls() )
library(dplyr)

source("defaults.R")
source("functions.R")
source("build_sim.R")

### simulate effects of scaling

args <- commandArgs(trailingOnly = TRUE)
j <- as.integer(args[1])

gxe_values <- c(0, 0.05)
#gxe_values <- c(0, 0.05, 0.1, 0.15)
alpha_values <- c(0.05,0.1,0.2,0.5)
LB_values <- seq(-1, 2, length = tuner)

# per-iteration output file
dir.create('Rdata/lam_box/', showWarnings=FALSE)
for(m in seq_along(gxe_values)){
  gxe_fixed <- gxe_values[m]
  for(n in seq_along(alpha_values)){
    alpha_fixed <- alpha_values[n]
    iter_savefile <- paste0('Rdata/lam_box/lam_box_gxe_', gxe_fixed,'_alpha_',alpha_fixed,'_iter_', j, '.Rdata')
    if(!file.exists(iter_savefile)){
      r2_iter <- array(NA, dim=c(length(LB_values),length(pgs_methods)),
                    dimnames=list(LB_values,pgs_methods))
      r2_diff_iter <- array(NA, dim=c(length(LB_values),length(pgs_methods[-1])),
                            dimnames=list(LB_values,pgs_methods[-1]))

      for(i in seq_along(LB_values)){
        temp <- run_sim(N, S, S_caus, h2_add, gxe_fixed,
                                   h2_coord, eta, alpha_fixed, LB_values[i])
        r2_iter[i,] <- temp[1:3]
        r2_diff_iter[i,] <- ((temp[2:3] - temp[1])/temp[1])*100
      }

      save( r2_iter, r2_diff_iter, j, LB_values, pgs_methods, file=iter_savefile)
    }
  }
}
