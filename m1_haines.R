# code from Haines 2019 Appendix 

rm(list = ls())

logit <- function(x){ log(x/(1-x)) }
expit <- function(x){ exp(x)/(1+exp(x)) }

##########################################
# Model 1: site covariates for abundance #
##########################################

# fit negative binomial

gM1nb <- function (param) {

  b0     <- (param[1])
  b1     <- (param[2])
  p      <- expit(param[3])
  r      <- exp(param[4])

  ptot <- 1 - (1 - p)^J

  xb <- b0 + b1 * xvec
  mu <- exp(xb)

  loglik <- 0

  for (i in 1:R) {

    term1 <- lgamma(r + yrow[i]) - lgamma(r)
    term2 <- r*log(r) + yrow[i] * xb[i]
    term3 <- yrow[i] * log(p) + sum(ymat[i, ] * jvec) * log(1 - p)
    term4 <- -(yrow[i] + r) * log(r + mu[i] * ptot)
    loglik <- loglik + term1 + term2 + term3 + term4

  }

  return(-loglik)
}

# fit poisson

gM1pc <- function (param) { 
  p       <- expit(param)
  ptot    <- 1 - (1 - p)^J
  loglikc <- ytot * log(p) + ysumj * log(1 - p) - ytot * log(ptot)

  return(-loglikc)
}


# generate with negbin 

gM1nbgen <- function (param) {

  b0 <- param[1]
  b1 <- param[2]
  p  <- param[3]
  r  <- param[4]
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



# setting

R <- 50
J <- 3

set.seed(125)

# set up

ymat <- matrix(0, R, J)
xvec <- read.table("xvec50.txt")
xvec <- as.matrix(xvec)

# parameters

pt  <- 0.615
b0t <- 3.077
b1t <- 0.101
rt  <- 11.415
rt  <- 4

# simulations

nsim   <- 10000
nprint <- 1000

b0vec  <- matrix(0, nsim, 1)
b1vec  <- matrix(0, nsim, 1)
pvec   <- matrix(0, nsim, 1)
pvecp  <- matrix(0, nsim, 1)
b0mvec <- matrix(0, nsim, 1)
b1mvec <- matrix(0, nsim, 1)
rvec   <- matrix(0, nsim, 1)

for (isim in 1:nsim) {

  # generate

  param <- c(b0t, b1t, pt, rt)
  ymat  <- gM1nbgen(param)

  # organize data

  yrow      <- rowSums(ymat)
  ycol      <- colSums(ymat)
  ytot      <- sum(sum(ymat))
  jvec      <- seq(0, J - 1)
  ysumj     <- sum(ycol * jvec)
  ylogfact  <- sum(sum(log(factorial(ymat))))
  yrlogfact <- sum(log(factorial(yrow)))

  # MLE

  param0      <- c(logit(pt))
  finopt      <- optim(param0, gM1pc, method="BFGS")
  pc          <- expit(finopt$par[1])
  pvecp[isim] <- pc

  param0      <- c(b0t, b1t, logit(pt), log(rt))
  finopt      <- optim(param0, gM1nb, method="Nelder-Mead")
  b0vec[isim] <- finopt$par[1]
  b1vec[isim] <- finopt$par[2]
  pvec[isim]  <- expit(finopt$par[3])
  rvec[isim]  <- exp(finopt$par[4])

  if (trunc(isim / nprint) * nprint == isim) {
    print(isim)
  }

}

#b0
mean(b0vec)
sd(b0vec)
quantile(b0vec,c(.025,.975))

#b1
mean(b1vec)
sd(b1vec)
quantile(b1vec,c(.025,.975))

#p
mean(pvec)
sd(pvec)
quantile(pvec,c(.025,.975))

#r
mean(rvec)
sd(rvec)
quantile(rvec,c(.025,.975))

#summaries
summary(b0vec)
summary(b1vec)
summary(pvec)
summary(rvec)

#histograms 
hist(b0vec)
hist(b1vec)
hist(pvec)
hist(rvec)


