rm( list=ls() )
library(dplyr)

source("defaults.R")
source("functions.R")
source("build_sim.R")

### simulate changes to coordinated GxC when uncoordinated GxC is decreased proportionally
args <- commandArgs(trailingOnly = TRUE)
j <- as.integer(args[1])

h2_coord_values <- seq(0, h2_gxe, length = tuner)

# per-iteration output file
dir.create('Rdata/h2_coord/', showWarnings=FALSE)
iter_savefile <- paste0('Rdata/h2_coord/h2_coord_iter_', j, '.Rdata')
if(!file.exists(iter_savefile)){
  r2_iter <- array(NA, dim=c(length(h2_coord_values),length(pgs_methods)),
                dimnames=list(h2_coord_values,pgs_methods))
  r2_diff_iter <- array(NA, dim=c(length(h2_coord_values),length(pgs_methods[-1])),
                    dimnames=list(h2_coord_values,pgs_methods[-1]))

  for(i in seq_along(h2_coord_values)){
    temp <- run_sim(N, S, S_caus, h2_add, h2_gxe - h2_coord_values[i],
                                   h2_coord_values[i], eta, alpha, lambda_boxcox)
    r2_iter[i,] <- temp[1:3]
    r2_diff_iter[i,] <- ((temp[2:3] - temp[1])/temp[1])*100
  }

  save( r2_iter, r2_diff_iter, j, h2_coord_values, pgs_methods, file=iter_savefile )
}


