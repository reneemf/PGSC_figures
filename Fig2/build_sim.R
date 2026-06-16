run_sim <- function(N, S, S_caus, h2_add, h2_uncoord, h2_coord, eta, alpha, lambda_boxcox, c_prob=0.5) {
  simdat  <- simulate_data(N*3, S, S_caus, h2_add, h2_uncoord, h2_coord, eta, alpha, lambda_boxcox, c_prob=c_prob)
  trn <- 1:N
  tst <- N+1:N
  val <- 2*N+1:N

  # run gwas
  gwas_out <- build_gwas(simdat$y[trn], simdat$G[trn,], simdat$C[trn])
  betas <- gwas_out["Estimate",]
  beta_pvs <- gwas_out[4,]

  # run gxewas
  gxewas_out <- build_gxewas(simdat$y[trn], simdat$G[trn,], simdat$C[trn])
  lambdas <- gxewas_out["Estimate",]
  lambda_pvs <- gxewas_out[4,]

  # run test pgs, ampPGS, PGSCs
  # find best pgs tau0
  test_pgs_R2s <- numeric(length(taus))
  names(test_pgs_R2s) <- taus

  for (i in seq_along(taus)){
    test_pgs <- build_pgs(simdat$G[tst,], betas, simdat$C[tst], beta_pvs, taus[i], simdat$y[tst])
    test_pgs_R2s[i] <- test_pgs$R2
  }
  pgs_tau_star <- taus[which.max(test_pgs_R2s)]
  opt_test_pgs <- build_pgs(simdat$G[tst,], betas, simdat$C[tst], beta_pvs, pgs_tau_star, simdat$y[tst])

  # build ampPGS w optimized tau0 to find rho
  test_ampPGS <- build_ampPGS(opt_test_pgs$PGS, betas, simdat$C[tst], beta_pvs, rho=NULL, pgs_tau_star, simdat$y[tst])

  # build PGSC w optimized tau0 to find rho(s)
  test_PGSC_R2s <- numeric(length(taus))
  for (i in seq_along(taus)){
    test_PGSC <- build_PGSC(opt_test_pgs$PGS, simdat$G[tst,], lambdas, simdat$C[tst], lambda_pvs,
                            rho=NULL, taus[i], simdat$y[tst])
    test_PGSC_R2s[i] <- test_PGSC$rho[["R2"]]
  }

  # grab best tau*
  PGSC_tau_star <- taus[which.max(test_PGSC_R2s)]

  # After finding tau star, re-run on test set at OPTIMAL tau to get correct rhos
  opt_test_PGSC <- build_PGSC(opt_test_pgs$PGS, simdat$G[tst,], lambdas, simdat$C[tst], lambda_pvs,
                              rho=NULL, PGSC_tau_star, simdat$y[tst])


  # run valid pgs, ampPGS, PGSCs
  valid_pgs <- build_pgs(simdat$G[val,], betas, simdat$C[val], beta_pvs, pgs_tau_star, simdat$y[val])
  valid_ampPGS <- build_ampPGS(valid_pgs$PGS, betas, simdat$C[val], beta_pvs,
                               test_ampPGS$rho, pgs_tau_star, simdat$y[val])
  valid_PGSC <- build_PGSC(valid_pgs$PGS, simdat$G[val,], lambdas, simdat$C[val], lambda_pvs,
                          opt_test_PGSC$rho, PGSC_tau_star, simdat$y[val])

  # calculate R2s
  pgs_R2 <- calc_R2s(simdat$y[val], valid_pgs$PGS,simdat$C[val])$R2
  ampPGS_R2 <- calc_R2s(simdat$y[val], valid_ampPGS$ampPGS,simdat$C[val])$R2
  PGSC_R2 <- calc_R2s(simdat$y[val], valid_PGSC$PGSC,simdat$C[val])$R2

  idx0 <- which(simdat$C[val] == min(simdat$C[val]))
  idx1 <- which(simdat$C[val] == max(simdat$C[val]))

  # Within-group R�
  pgs0_R2  <- summary(lm(simdat$y[val][idx0] ~ valid_pgs$PGS[idx0]))$r.squared
  pgs1_R2  <- summary(lm(simdat$y[val][idx1] ~ valid_pgs$PGS[idx1]))$r.squared
  PGSC0_R2 <- summary(lm(simdat$y[val][idx0] ~ valid_PGSC$PGSC[idx0]))$r.squared
  PGSC1_R2 <- summary(lm(simdat$y[val][idx1] ~ valid_PGSC$PGSC[idx1]))$r.squared

  c( pgs_R2, ampPGS_R2, PGSC_R2,
    pgs0_R2, pgs1_R2, PGSC0_R2, PGSC1_R2 )

}


