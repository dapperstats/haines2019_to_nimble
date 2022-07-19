# code from Haines 2019 Appendix 

##########################################
# Model 1: site covariates for abundance #
##########################################

# fit negative binomial

gM1nb <- function (param, ymat, J, R, xvec) {

  yrow     <- as.matrix(rowSums(ymat))
  ycol     <- as.matrix(colSums(ymat))
  ytot     <- sum(sum(ymat))
  jvec     <- seq(0, J - 1)
  ysumj    <- sum(ycol * jvec)
  ylogfact <- sum(sum(log(factorial(ymat))))



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

gM1nbgen <- function (param, J, R, xvec) {

  b0 <- param[1]
  b1 <- param[2]
  p  <- expit(param[3])
  r  <- exp(param[4])

  ymat <- matrix(0, R, J + 1)
  prob <- numeric(J + 1)
  for (i in 1:(J)) {
    prob[i] <- pow(1 - p, i - 1) * p
  }
  prob[J + 1] <- 1 - sum(prob[1:J])

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


