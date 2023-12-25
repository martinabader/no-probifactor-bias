
# ---------------------------------------------------------------------------------------------
#
# title: "No probifactor model fit index bias"
# analysis: "replication simulation population level"
#
# ---------------------------------------------------------------------------------------------



library(lavaan)


# functions ---------------------------------------------------------------


# minimum of ML fit function 
calcFmin <- function(SigmaHat, S){
  fmin <- sum(diag(S %*% solve(SigmaHat))) + log(det(SigmaHat)) - log(det(S)) - ncol(S)
  return(fmin)
} 


# population RMSEA
getRMSEA.pop <- function(Fmin, df){
  rmsea <- sqrt(Fmin/df)
  return(rmsea)
}



# basic population model --------------------------------------------------


## note:
# loadings of all indicators - except for indicator y4 - are in standardized metric


#  loadings
lambda <- matrix(ncol=3, byrow = T, c(
  c(0.850,0.000,0.000), 
  c(0.860,0.000,0.000), 
  c(0.850,0.000,0.000),
  c(0.000,0.810,0.000),  # in cross-loading conditions, this is no longer std, b/c residual was always fixed to 1- sum(lambda[4,]^2) by Greene et al. 
  c(0.000,0.600,0.000),
  c(0.000,0.740,0.000),
  c(0.000,0.000,0.670),
  c(0.000,0.000,0.660),
  c(0.000,0.000,0.790),
  c(0.000,0.000,0.880),
  c(0.000,0.000,0.870)
))

# latent correlations
phi <- diag(1,ncol(lambda))
phi[1,2] <- phi[2,1] <- .82
phi[2,3] <- phi[3,2] <- .60
phi[1,3] <- phi[3,1] <- .59


# residuals
theta <- diag(diag(1 - lambda %*% phi %*% t(lambda)))




# conditions --------------------------------------------------------------

# 1. no model complexities
# 2. ustd. cross-factor loading of .10
# 3. ustd. cross-factor loading of .30
# 4. ustd. cross-factor loading of .50 
# 5. correlated residual (within factor) of .10
# 6. correlated residual (within factor) of .30
# 7. correlated residual (within factor) of .50
# 8. correlated residual (between factors) of .10
# 9. correlated residual (between factors) of .30
# 10. correlated residual (between factors) of .50


# number of observations for BIC
N <- 10000




# simulation misspecified models ------------------------------------------

# results object

colnames.res <- c('condition', 
                  'fmin.cf', 'df.cf', 'npar.cf', 'rmsea.cf', 'bic.cf',
                  'fmin.ho', 'df.ho', 'npar.ho', 'rmsea.ho', 'bic.ho',
                  'fmin.bf', 'df.bf', 'npar.bf', 'rmsea.bf', 'bic.bf')

res <- as.data.frame(matrix(NA, ncol = length(colnames.res), nrow = 10))

colnames(res) <- colnames.res




for(i in 1:10){
  
  
  ### define condition ###
  
  c.cond <- i
  
  
  ## notes:
  # residual correlations are defined in unstandardized metric
  # cross-loadings are defined in unstandardized metric
  
  if(c.cond == 1){cross.l <- 0; cor.res.w <- 0; cor.res.b <- 0
  } else if(c.cond == 2){cross.l <- .10; cor.res.w <- 0; cor.res.b <- 0             # std. cross-factor loading of .094
  } else if(c.cond == 3){cross.l <- .30; cor.res.w <- 0; cor.res.b <- 0             # std. cross-factor loading of .254
  } else if(c.cond == 4){cross.l <- .50; cor.res.w <- 0; cor.res.b <- 0             # std. cross-factor loading of .388
  } else if(c.cond == 5){cross.l <- .00; cor.res.w <- .05577104; cor.res.b <- 0     # std. within-factor correlated residual of .10
  } else if(c.cond == 6){cross.l <- .00; cor.res.w <- .1673131; cor.res.b <- 0      # std. within-factor correlated residual of .30
  } else if(c.cond == 7){cross.l <- .00; cor.res.w <- .2788552; cor.res.b <- 0      # std. within-factor correlated residual of .50
  } else if(c.cond == 8){cross.l <- .00; cor.res.w <- .0; cor.res.b <- .03089211    # std. between-factor correlated residual of .10
  } else if(c.cond == 9){cross.l <- .00; cor.res.w <- .0; cor.res.b <- .09267633    # std. between-factor correlated residual of .30
  } else if(c.cond == 10){cross.l <- .00; cor.res.w <- .0; cor.res.b <- .1544606}   # std. between-factor correlated residual of .50
  
  
  # cross-loading?
  lambda[4, 1] <- cross.l
  
  # correlated residual within factor?
  theta[7,8] <- theta[8,7] <- cor.res.w
  
  # correlated residual between factors?
  theta[3,4] <- theta[4,3] <- cor.res.b
  
  # replicate Greene's approach of defining the residual of y4 
  theta[4,4] <- 1- sum(lambda[4,]^2)
  
  
  

  ### create population covariance matrix ###
  
  sigma <- lambda %*% phi %*% t(lambda) + theta
  rownames(sigma) <- colnames(sigma) <- paste0('y', 1:nrow(lambda))
  
  
  
  
  ### analyze data with different analysis models ###
  
  ## pure correlated factors model
  
  m.cf <- '
  f1 =~ NA*y1 + y2 + y3 
  f2 =~ NA*y4 + y5 + y6
  f3 =~ NA*y7 + y8 + y9 + y10 + y11
  f1 ~~ 1*f1
  f2 ~~ 1*f2
  f3 ~~ 1*f3
  '
  res.cf <- sem(m.cf, sample.cov = sigma, sample.nobs = N, sample.cov.rescale=F) 

  # get minimum of fit function
  cov.cf <- fitted(res.cf)$cov # model-implied covariance matrix
  Fmin.cf <- calcFmin(cov.cf, sigma)
  names(Fmin.cf) <- 'Fmin'
  
  # get population RMSEA
  df.cf <-  fitmeasures(res.cf)['df']
  RMSEA.cf <- getRMSEA.pop(Fmin.cf, df.cf)
  
  # get number of free parameters
  k.cf <- fitmeasures(res.cf)['npar']

  # get BIC
  BIC.cf <- fitmeasures(res.cf)['bic']
  
  fit.cf <- c(Fmin.cf, df.cf, k.cf, RMSEA.cf, BIC.cf)
              
              
  
  
  ## pure higher-order model
  
  m.ho <- '
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
  f3 =~ y7 + y8 + y9 + y10 + y11
  
  g =~ NA*f1 + f2 + f3
  
  g ~~ 1*g
  '
  res.ho <- sem(m.ho, sample.cov = sigma, sample.nobs = N, sample.cov.rescale=F)

  # get minimum of fit function
  cov.ho <- fitted(res.ho)$cov # model-implied covariance matrix
  Fmin.ho <- calcFmin(cov.ho, sigma)
  names(Fmin.ho) <- 'Fmin'

  # get population RMSEA
  df.ho <-  fitmeasures(res.ho)['df']
  RMSEA.ho <- getRMSEA.pop(Fmin.ho, df.ho)
  
  # get number of free parameters
  k.ho <- fitmeasures(res.ho)['npar']
  
  # get BIC
  BIC.ho <- fitmeasures(res.ho)['bic']
  
  fit.ho <- c(Fmin.ho, df.ho, k.ho, RMSEA.ho, BIC.ho)
  
  
  
  
  ## pure bifactor model
  
  m.bf <- '
  g =~ NA*y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10 + y11
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
  f3 =~ y7 + y8 + y9 + y10 + y11
  
  g ~~ 1*g
  
  f1 + f2 + f3 ~~ 0*g
  f2 + f3 ~~ 0*f1
  f3 ~~ 0*f2
  '
  
  res.bf <- sem(m.bf, sample.cov = sigma, sample.nobs = N, sample.cov.rescale=F)

  # get minimum of fit function
  cov.bf <- fitted(res.bf)$cov # model-implied covariance matrix
  Fmin.bf <- calcFmin(cov.bf, sigma)
  names(Fmin.bf) <- 'Fmin'

  # get population RMSEA
  df.bf <-  fitmeasures(res.bf)['df']
  RMSEA.bf <- getRMSEA.pop(Fmin.bf, df.bf)
  
  # get number of free parameters
  k.bf <- fitmeasures(res.bf)['npar']
  
  # get BIC 
  BIC.bf <- fitmeasures(res.bf)['bic']
  
  fit.bf <- c(Fmin.bf, df.bf, k.bf, RMSEA.bf, BIC.bf)
  
  
  
  ### save results ###
  
  res[i, ] <- c(i, fit.cf, fit.ho, fit.bf)
  
  
  
  ## go to next condition
  
  i <- i + 1

}




round(res, 4)


















