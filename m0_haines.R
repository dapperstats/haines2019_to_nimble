# code from Haines 2019 Appendix 

rm(list = ls())

logit <- function(x){ log(x/(1-x)) }
expit <- function(x){ exp(x)/(1+exp(x)) }


####################################
# Model 0: constant mu, constant p #
####################################
# note that the "yrlogfact" has been replaced with "ylogfact" to match 
# the math and mathematica code

# fit negative binomial
 
gM0nb <- function (param) {

  mu     <- exp(param[1])
  p      <- expit(param[2])
  r      <- exp(param[3])

  ptot   <- 1 - (1 - p)^J

  term1  <- sum(lgamma(r + yrow)) - R * lgamma(r) - ylogfact
  term2  <- R * r * log(r) + ytot * log(mu)
  term3  <- ytot * log(p) + ysumj * log(1-p)
  term4  <- -(ytot + R * r) * log(r + mu * ptot)
  loglik <- term1 + term2 + term3 + term4

  return(-loglik)

}


# generate with negbin mixing
# edited to transform the parameters so as to align with density function
gM0nbgen <- function (param) {

  mu   <- exp(param[1])
  p    <- expit(param[2])
  r    <- exp(param[3])
  ymat <- matrix(0, R, J+1)

  prob <- c(p, (1 - p) * p, (1 - p)^2 * p, (1 - p)^3)

  for (i in 1:R) {
    n <- rnbinom(1, size = r, mu = mu)
    if (n > 0) {
      ymat[i, ] <- rmultinom(1, n, prob)
    }
  }
  ymat <- ymat[ , 1:J]

  return(ymat)
}


