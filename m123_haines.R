
# combining

##########################################
# Model 1: site covariates for abundance #
##########################################

# and

##########################################
# Model 2: site covariates for capture   #
##########################################

# and

#############################################
# Model 3: removal covariates for capture   #
#############################################


gM123nb <- function (param, ymat, J, R, xvec, tvecR, tvecJ) {

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
  g2 <- param[5]
  r  <- exp(param[6])

  mu <- exp(b0 + b1 * xvec)

  loglik <- -ylogfact

  for (i in 1:R) {

    p <- expit(g0 + g1 * tvecR[i] + g2 * tvecJ)

    pp <- c(0, p)


    prob <- numeric(J)
    for (j in 1:J) {
      prob[j] <- prod(1 - pp[1:j]) * p[j]

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

gM123nbgen <- function (param, J, R, xvec, tvecR, tvecJ) {

  b0 <- param[1]
  b1 <- param[2]
  g0 <- param[3]
  g1 <- param[4]
  g2 <- param[5]
  r  <- exp(param[6])

  mu <- exp(b0 + b1 * xvec)

  ymat <- matrix(0, R, J + 1)

  for(i in 1:R) {
    
    n   <- rnbinom(1, size = r, mu = mu[i])
    if(n > 0) {
      p <- expit(g0 + g1 * tvecR[i] + g2 * tvecJ)

      pp <- c(0, p)
      prob <- numeric(J + 1)
      for (j in 1:J) {

        prob[j] <-  prod(1 - pp[1:j]) * p[j]

      }
      prob[J + 1] <- 1 - sum(prob[1:J])

      ymat[i,] <- rmultinom(1, n, prob)
    }
  }
  ymat <- ymat[,1:J]

  return(ymat)
}


