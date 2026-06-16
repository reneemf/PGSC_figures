simulate_data <- function(N,S,S_caus,h2_add,h2_uncoord,h2_coord,eta,alpha,lambda_boxcox,c_prob=0.5){
  # genotype & context
  G	<- scale( replicate( S, rbinom( N, 2, runif(1,0.05,0.5) ) ) ) # MAF 5-50%
  caus  <- sample(S,S_caus,rep=F)
  C_raw <- rbinom(N,1,c_prob)
  C <- scale(C_raw)

  # construct gamma (h2_coord = gamma^2 * h2_add)
  gamma <- sqrt(h2_coord/h2_add)

  # std devs # need to add Cs
  sigma_epsilon <- sqrt(1 - h2_add - h2_coord - h2_uncoord)
  sigma_f <- sqrt((2*sigma_epsilon^2)/(1+eta))
  sigma_m <- sqrt(eta*sigma_f^2)

  # effects
  beta   <- rep(0,S)
  lambda <- rep(0,S)
  beta[caus]   <- sqrt(h2_add/S_caus)    *scale(rnorm(S_caus))
  lambda[caus] <- sqrt(h2_uncoord/S_caus)*scale(rnorm(S_caus))

  epsilon <- numeric(N)
  epsilon[C_raw == 0] <- rnorm(sum(C_raw == 0), 0, sigma_f)
  epsilon[C_raw == 1] <- rnorm(sum(C_raw == 1), 0, sigma_m)

  # phenotype modeling
  y <- G%*%beta + gamma*C*(G%*%beta) + C*(G%*%lambda) + C*sqrt(alpha) + epsilon

  # Lambda Box-Cox: shift y to be >= 1, then transform. Handles lambda=0 via log.
  if(lambda_boxcox != 1){
    y <- y - min(y) + 1
    if (lambda_boxcox == 0){
      y <- log(y)
    } else {
      y <- (y^lambda_boxcox - 1) / lambda_boxcox
    }
    y <- scale(y)
  }

  list(y=y, G=G, C=C)
}

### build the gwas
# input: phenotypes, training genotypes, context
build_gwas <- function(y, G, C){
  apply( G, 2, function(G) summary( lm( y ~ G + C) )$coef['G',] )
}

### build the gxewas
build_gxewas <- function(y, G, C){
  apply( G, 2, function(G) summary( lm( y ~ G + C + G:C) )$coef['G:C',] )
}

### build the pgs
# input: testing/validation genotypes, betas, pvs, pv thresh
build_pgs <- function(G, B, C, pv, thresh, y){
  B[pv > thresh] <- 0
  PGS <- G%*%B
  lm_out <- calc_R2s(y,PGS,C)
  list( PGS=PGS, R2=lm_out$R2 )
}

### build the ampPGS
build_ampPGS <- function(pgs, B, C, pv, rho, thresh, y){
  B[pv > thresh] <- 0
  if( is.null(rho) ){
    rho <- learn_rho( y, C, pgs, pgs )
  }
  ampPGS <- pgs + rho$rho * pgs * C
  list( ampPGS=ampPGS, B=B, rho=rho)
}

### build the PGSC
build_PGSC <- function(pgs, G, L, C, pv, rho, thresh, y){
  L[pv > thresh] <- 0
  PGxCS <- G%*%L
  PGxCS_min <- PGxCS*as.numeric(C==min(C))
  PGxCS_max <- PGxCS*as.numeric(C==max(C))
  if(is.null(rho)){
    rho <- learn_rho( y, C, pgs, PGxCS_min, PGxCS_max )
  }
  PGSC <- pgs + rho$rho1 * PGxCS_min + rho$rho2 * PGxCS_max
  list(PGSC=PGSC, rho=rho, PGxCS=PGxCS, PGxCS_min=PGxCS_min, PGxCS_max=PGxCS_max)
}

### learn tuning params
learn_rho <- function( y, C, x1, x2, x3=NULL ){
  if(is.null(x3)){
    model_out <- lm( y ~ x1 + x2:C + C )
    weights <- coef( model_out )[c('x1','x2:C')]
    R2 <- summary(model_out)$r.squared
    if (is.na(weights[1]) || abs(weights[1]) < 1e-10) {
      return(list(rho = 0, R2 = R2))
    }
    rho <- as.numeric(weights[2]/weights[1])
    return(list(rho=as.numeric(weights[2]/weights[1]), R2=R2))
  }else{
    model_out <- lm( y ~ x1 + x2 + x3 + C )
    weights <- coef( model_out )[c('x1','x2','x3')]
    R2 <- summary(model_out)$r.squared
    if (is.na(weights[1]) || abs(weights[1]) < 1e-10) {
      return(list(rho1 = 0, rho2 = 0, R2 = R2))
    }
    return(list(rho1=as.numeric(weights[2]/weights[1]),
                rho2=as.numeric(weights[3]/weights[1]), R2=R2))
  }
}

### calculate the R2
calc_R2s <- function(y, pgs,C){
  y <- resid( lm( y ~ C ) )
  pgs <- resid( lm( pgs ~ C ) )
  lm_out <- summary(lm(y ~ pgs))
  R2 <- lm_out$r.squared
  list(lm_out=lm_out, R2=lm_out$r.squared)
}



