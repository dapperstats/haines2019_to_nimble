# code from Haines 2019 Appendix 

rm(list = ls())

logit <- function(x){ log(x/(1-x)) }
expit <- function(x){ exp(x)/(1+exp(x)) }


####################################
# Model 0: constant mu, constant p #
####################################

# fit negative binomial
 
gM0nb <- function (param) {

  mu     <- exp(param[1])
  p      <- expit(param[2])
  r      <- exp(param[3])

  ptot   <- 1 - (1 - p)^J

  term1  <- sum(lgamma(r + yrow)) - R * lgamma(r) - yrlogfact
  term2  <- R * r * log(r) + ytot * log(mu)
  term3  <- ytot * log(p) + ysumj * log(1-p)
  term4  <- -(ytot + R * r) * log(r + mu * ptot)
  loglik <- term1 + term2 + term3 + term4

  return(-loglik)

}

# generate with negbin mixing

gM0nbgen <- function (param) {

  mu <- param[1]
  p  <- param[2]
  r  <- param[3]
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



# negative binomial bootstraps
# setting

#R  <- 20
#J  <- 3

#set.seed <- 125

# parameters

#pt  <- 0.615
#mut <- 28.048
#rt  <- 11.415

# simulations

#nsim   <- 10000
#nprint <- 1000
#muvec  <- matrix(0,nsim,1)
#pvec   <- matrix(0,nsim,1)
#rvec   <- matrix(0,nsim,1)

#for (isim in 1:nsim) {

  # generate

#  param <- c(mut, pt, rt)
#  ymat  <- gM0nbgen(param)

  # organize data

#  yrow      <- as.matrix(rowSums(ymat))
#  ycol      <- as.matrix(colSums(ymat))
#  ytot      <- sum(sum(ymat))
#  jvec      <- seq(0, J - 1)
#  ysumj     <- sum(ycol * jvec)
#  ylogfact  <- sum(sum(log(factorial(ymat))))
#  yrlogfact <- sum(log(factorial(yrow)))

  # MLE

#  param0 <- c(log(mut), logit(pt), log(rt))
#
#  finopt <- optim(param0, gM0nb, method="Nelder-Mead")
#
#  muvec[isim] <- exp(finopt$par[1])
#  pvec[isim]  <- expit(finopt$par[2])
#  rvec[isim]  <- exp(finopt$par[3])
#
#  if (trunc(isim / nprint) * nprint == isim) {
#    print(isim)
#  }
#
#}




