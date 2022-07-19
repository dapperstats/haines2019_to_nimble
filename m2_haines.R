# code from Haines 2019 Appendix 

##########################################
# Model 2: site covariates for capture   #
##########################################
# note that the ylogfact is moved to pre-loop, as in Model 3, 
# the term3 calculation is represented to be similar to the other models
# term3 <- sum(ymat[i,] * log(prob)) # identical calculation (this is what's in the ms)

# fit negative binomial

gM2nb <- function (param, ymat, J, R, tvec) {

  yrow     <- as.matrix(rowSums(ymat))
  ycol     <- as.matrix(colSums(ymat))
  ytot     <- sum(sum(ymat))
  jvec     <- seq(0, J - 1)
  ysumj    <- sum(ycol * jvec)
  ylogfact <- sum(sum(log(factorial(ymat))))


  mu   <- exp(param[1])
  r    <- exp(param[4])

  g0 <- param[2]
  g1 <- param[3]
  p  <- expit(g0 + g1 * tvec)



  loglik <- -ylogfact

  for (i in 1:R) {
 
    ptot <- 1 - (1 - p[i])^J
   
    prob <- numeric(J)
    for (j in 1:J) {
      prob[j] <- pow(1 - p[i], j - 1) * p[i]
    }

    term1 <- lgamma(r + yrow[i]) - lgamma(r)
    term2 <- r * log(r) + yrow[i] * log(mu)
    term3 <- yrow[i] * log(p[i]) + sum(ymat[i, ] * jvec) * log(1 - p[i])
    term4 <- -(yrow[i] + r) * log(r + mu * ptot)
    loglik <- loglik + term1 + term2 + term3 + term4

  }

  return(-loglik)
}



# generate with negbin 
# edited to transform the parameters so as to align with density function


gM2nbgen <- function (param, J, R, tvec) {

  mu   <- exp(param[1])
  r    <- exp(param[4])

  g0 <- param[2]
  g1 <- param[3]
  p  <- expit(g0 + g1 * tvec)

  ymat <- matrix(0, R, J + 1)

  for(i in 1:R) {
    n   <- rnbinom(1, size = r, mu = mu)
    if(n > 0) {
      probi <- numeric(J + 1)
      for (j in 1:J) {
        probi[j] <- pow(1 - p[i], j - 1) * p[i]
      }
      probi[J + 1] <- 1 - sum(probi[1:J])
      ymat[i,] <- rmultinom(1, n, probi)
    }
  }
  ymat <- ymat[,1:J]

  return(ymat)
}



