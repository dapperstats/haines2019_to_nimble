# code from Haines 2019 Appendix 

rm(list = ls())

logit <- function(x){ log(x/(1-x)) }
expit <- function(x){ exp(x)/(1+exp(x)) }

##########################################
# Model 1: site covariates for capture   #
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

gM2nbgen <- function (param) {

  b0 <- param[1]
  b1 <- param[2]
  p  <- expit(param[3])
  r  <- exp(param[4])
  ymat <- matrix(0, R, J + 1)
  prob <- c(p, (1 - p) * p,(1 - p)^2 * p, (1 - p)^3)

  for(i in 1:R) {
    mui <- exp(b0 + b1 * xvec[i])
    n   <- rnbinom(1, size = r, mu = mui)
    if(n > 0) {
      ymat[i,] <- rmultinom(1, n, prob)
    }
  }
  ymat <- ymat[,1:J]

  return(ymat)
}



