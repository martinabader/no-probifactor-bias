
# ---------------------------------------------------------------------------------------------
#
# title: "No probifactor model fit index bias"
# analysis: "Nesting and equivalence testing (Bentler & Satorra, 2010)"
#
# ---------------------------------------------------------------------------------------------


library(lavaan)



# functions ---------------------------------------------------------------

# generates normally distributed data given cov and N using Cholesky decomposition

genData <- function(SigmaHat, N){
  p <- ncol(SigmaHat)
  randomData <- matrix(rnorm(N*p), N, p) 
  
  cc <- chol(SigmaHat)
  out <- t(t(cc) %*% t(randomData))
  return(out)
}


# get minimum of ML fitting function 

calcFmin <- function(SigmaHat, S){
  fmin <- sum(diag(S %*% solve(SigmaHat))) + log(det(SigmaHat)) - log(det(S)) - ncol(S)
  return(fmin)
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

conditions <- seq(1,10)




# container ---------------------------------------------------------------

# object to save results

colnames.res <- c('condition', 'fmin.ho', 'fmin.bf')

res <- as.data.frame(matrix(NA, ncol = length(colnames.res), nrow = 10))

colnames(res) <- colnames.res




# NET ---------------------------------------------------------------------


for(i in 1:10){
  
  
  ## define condition
  
  c.cond <- conditions[i]
  
  if(c.cond == 1){cross.l <- 0; cor.res.w <- 0; cor.res.a <- 0
  } else if(c.cond == 2){cross.l <- .10; cor.res.w <- 0; cor.res.a <- 0
  } else if(c.cond == 3){cross.l <- .30; cor.res.w <- 0; cor.res.a <- 0
  } else if(c.cond == 4){cross.l <- .50; cor.res.w <- 0; cor.res.a <- 0
  } else if(c.cond == 5){cross.l <- .00; cor.res.w <- .05577104; cor.res.a <- 0
  } else if(c.cond == 6){cross.l <- .00; cor.res.w <- .1673131; cor.res.a <- 0  
  } else if(c.cond == 7){cross.l <- .00; cor.res.w <- .2788552; cor.res.a <- 0 
  } else if(c.cond == 8){cross.l <- .00; cor.res.w <- .0; cor.res.a <- .03089211
  } else if(c.cond == 9){cross.l <- .00; cor.res.w <- .0; cor.res.a <- .09267633
  } else if(c.cond == 10){cross.l <- .00; cor.res.w <- .0; cor.res.a <- .1544606}
  
  # cross-loading?
  lambda[4, 1] <- cross.l
  
  # correlated residual within factor?
  theta[7,8] <- theta[8,7] <- cor.res.w
  
  # correlated residual across factors?
  theta[3,4] <- theta[4,3] <- cor.res.a
  
  # replicate std error for residual of y4
  theta[4,4] <- 1- sum(lambda[4,]^2)
  
  
  
  ## get population sigma
  
  sigma <- lambda %*% phi %*% t(lambda) + theta
  rownames(sigma) <- colnames(sigma) <- paste0('y', 1:nrow(lambda))
  
  
  ## generate data 
  
  set.seed(1111) # seed for replicability
  c.data <- genData(sigma, N = 1000)
  
  
  ## analyses models ##
  
  # if a complexitiy (such as a cross-loading) is present in the population model,
  # the analyses models include the respective model complexity as well
  
  if(c.cond == 1){c.complexity <- NULL  # no model complexity
  } else if(c.cond %in% c(2,3,4)){c.complexity <- 'f1 =~ NA*y4'  # cross-loading
  } else if(c.cond %in% c(5,6,7)){c.complexity <- 'y7 ~~ y8'     # within-factor correlated residual
  } else if(c.cond %in% c(8,9,10)){c.complexity <- 'y3 ~~ y4'}   # between-factor correlated residual
  
  

    ## correlated factors model
    m.cf <- '
    f1 =~ NA*y1 + y2 + y3
    f2 =~ NA*y4 + y5 + y6
    f3 =~ NA*y7 + y8 + y9 + y10 + y11
    
    f1 ~~ 1*f1
    f2 ~~ 1*f2
    f3 ~~ 1*f3
    '
    
    # add complexity of current condition to model syntax
    m.cf <- paste0(m.cf, c.complexity) 
    
    # fit model to sample data
    res.cf <- sem(m.cf, data = c.data)
    
    # get model-implied covariance matrix
    cov.cf <- fitted(res.cf)$cov
    

    
    
    ## higher-order model
    m.ho <- '
    f1 =~ y1 + y2 + y3
    f2 =~ y4 + y5 + y6
    f3 =~ y7 + y8 + y9 + y10 + y11
    
    g =~ NA*f1 + f2 + f3
    
    g ~~ 1*g
    '
    
    # add complexity of current condition to model syntax
    m.ho <- paste0(m.ho, c.complexity) 
    
    # fit model (based on covariance matrix implied by the correlated factors model)
    res.ho <- sem(m.ho, sample.cov = cov.cf, sample.nobs = 1000, sample.cov.rescale=F)
    
    # get model-implied covariance matrix
    cov.ho <- fitted(res.ho)$cov
    
    # get minimum of fit function
    Fmin.ho <- calcFmin(cov.ho, cov.cf)

    
    
    
    ## bifactor model
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
    
    # add complexity of current condition to model syntax
    m.bf <- paste0(m.bf, c.complexity) 
    
    # fit model (based on covariance matrix implied by the correlated factors model)
    res.bf <- sem(m.bf, sample.cov = cov.cf, sample.nobs = 1000, sample.cov.rescale=F)
    
    # get model-implied covariance matrix
    cov.bf <- fitted(res.bf)$cov
    
    # get minimum of fit function
    Fmin.bf <- calcFmin(cov.bf, cov.cf)

    
    
    ### save results ###
    
    res[i, ] <- c(i, Fmin.ho, Fmin.bf)
    
    
    ## go to next condition
    
    i <- i + 1
 

}


# Guideline by Asparouhov & Muthén (2019, p. 303; doi: 10.1080/10705511.2018.1513795):
# "Because F0 is computed numerically through the minimization
# procedure, precise zero values are unrealistic. Instead,
# small values are used as cutoff values and 10^-7 appears to
# work well enough for most situations."


res

# => covariance equivalence
















