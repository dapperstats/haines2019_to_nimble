# code from Haines 2019 Appendix 

rm(list = ls())

logit <- function(x){ log(x/(1-x)) }
expit <- function(x){ exp(x)/(1+exp(x)) }

##########################################
# Model 2: site covariates for capture   #
##########################################

# fit negative binomial

gM2nb <- function (param) {

  b0     <- (param[1])
  b1     <- (param[2])
  p      <- expit(param[3])
  r      <- exp(param[4])

  ptot <- 1 - (1 - p)^J

  xb <- b0 + b1 * xvec
  mu <- exp(xb)

  loglik <- -ylogfact

  for (i in 1:R) {

    term1 <- lgamma(r + yrow[i]) - lgamma(r)
    term2 <- r*log(r) + yrow[i] * xb[i]
    term3 <- yrow[i] * log(p) + sum(ymat[i, ] * jvec) * log(1 - p)
    term4 <- -(yrow[i] + r) * log(r + mu[i] * ptot)
    loglik <- loglik + term1 + term2 + term3 + term4

  }

  return(-loglik)
}



# generate with negbin 
# edited to transform the parameters so as to align with density function

gM2nbgen <- function (param) {

  mu   <- exp(param[1])
  r    <- exp(param[4])

  g0t <- param[2]
  g1t <- param[3]
  p  <- expit(g0t + g1t * tvec)

  ymat <- matrix(0, R, J + 1)

  for(i in 1:R) {
    n   <- rnbinom(1, size = r, mu = mu)
    if(n > 0) {
      probi <- c(p[i], (1 - p[i]) * p[i], (1 - p[i])^2 * p[i], (1 - p[i])^3)

      ymat[i,] <- rmultinom(1, n, probi)
    }
  }
  ymat <- ymat[,1:J]

  return(ymat)
}



