
# combining

##########################################
# Model 1: site covariates for abundance #
##########################################

# and

##########################################
# Model 2: site covariates for capture   #
##########################################




gM12nb <- function (param, ymat, J, R, xvec, tvec) {

  yrow     <- as.matrix(rowSums(ymat))
  ycol     <- as.matrix(colSums(ymat))
  ytot     <- sum(sum(ymat))
  jvec     <- seq(0, J - 1)
  ysumj    <- sum(ycol * jvec)
  ylogfact <- sum(sum(log(factorial(ymat))))



  b0 <- (param[1])
  b1 <- (param[2])
  g0 <- param[3]
  g1 <- param[4]
  r  <- exp(param[5])

  p  <- expit(g0 + g1 * tvec)
  mu <- exp(b0 + b1 * xvec)

  loglik <- -ylogfact

  for (i in 1:R) {


    prob <- numeric(J)
    for (j in 1:J) {
      prob[j] <- pow(1 - p[i], j - 1) * p[i]
    }
    ptot <- sum(prob)

    term1 <- lgamma(r + yrow[i]) - lgamma(r)
    term2 <- r*log(r) + yrow[i] * log(mu[i])
    term3 <- sum(ymat[i,] * log(prob))
    term4 <- -(yrow[i] + r) * log(r + mu[i] * ptot)
    loglik <- loglik + term1 + term2 + term3 + term4

  }

  return(-loglik)
}



# generate with negbin 

gM12nbgen <- function (param, J, R, xvec, tvec) {

  b0 <- param[1]
  b1 <- param[2]
  g0 <- param[3]
  g1 <- param[4]

  r  <- exp(param[5])
  p  <- expit(g0 + g1 * tvec)
  mu <- exp(b0 + b1 * xvec)

  ymat <- matrix(0, R, J + 1)

  for(i in 1:R) {
    
    n   <- rnbinom(1, size = r, mu = mu[i])
    if(n > 0) {

      prob <- numeric(J + 1)
      for (j in 1:J) {
        prob[j] <- pow(1 - p[i], j - 1) * p[i]
      }
      prob[J + 1] <- 1 - sum(prob[1:J])

      ymat[i,] <- rmultinom(1, n, prob)
    }
  }
  ymat <- ymat[,1:J]

  return(ymat)
}


