
# multinomial negative binomial mixture model
#   assuming a rectangular  matrix of observations
#   but generalizing to allow more flexibles


dM3_nb <- nimbleFunction (
  run = function (x   = double(2),
                  mut = double(),
                  g0  = double(),
                  g1  = double(),
                  rt  = double(),
                  J   = integer(),
                  R   = integer(),
                  w   = double(1), 
                  log = logical(0, default = 0)) {

    mu   <- exp(mut)
    pt   <- g0 + g1 * w
    p    <- expit(pt)
    r    <- exp(rt)
  pp <- c(0, p)

  x_row      <- integer(R)
  x_miss_row <- integer(R)

  for (i in 1:R) {
    x_row[i]      <- sum(x[i, ])
    x_miss_row[i] <- x[i, ] %*% seq(0, J - 1)
  }

  x_tot  <- sum(x)
  x_vec  <- seq(0, J - 1)
  x_sumj <- sum(x_miss_row[1:R]) 

  x_logfact     <- sum(lfactorial(x))
  
  logProb <- -x_logfact

  for (i in 1:R) {
 
  
    prob <- double(J)
    for (j in 1:J) {
      prob[j] <- prod(1 - pp[1:j]) * p[j]
    }
    ptot <- sum(prob)

    term1   <- lgamma(r + x_row[i]) - lgamma(r)
    term2   <- r * log(r) + x_row[i] * log(mu)
    term3   <- sum(x[i, ] * log(prob))
    term4   <- -(x_row[i] + r) * log(r + mu * ptot)
    logProb <- logProb + term1 + term2 + term3 + term4

  }
  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())  
})


rM3_nb <- nimbleFunction(
  run = function(n   = integer(),
                 mut = double(),
                 g0  = double(),
                 g1  = double(),
                 rt  = double(),
                 J   = integer(), 
                 R   = integer(),
                 w   = double(1)) {

    mu   <- exp(mut)
    pt   <- g0 + g1 * w
    p    <- expit(pt)
    r    <- exp(rt)
    
    ans <- matrix(0, nrow = R, ncol = J + 1)

    pp <- c(0, p)

    for (i in 1:R) {
      n <- rnbinom(n = 1, size = r, mu = mu)
      if (n > 0) {
        prob <- double(J + 1)
        for (j in 1:J) {
          prob[j] <- prod(1 - pp[1:j]) * p[j]
        }
        prob[J + 1] <- 1 - sum(prob[1:J])


        ans[i, ] <- rmulti(n = 1, size = n, prob = prob)
      }
    }

    return(ans[ , 1:J])
    returnType(integer(2))
})

registerDistributions(list(
  dM3_nb = list(
    BUGSdist = "dM3_nb(mut, g0, g1, rt, J, R, z)",
    Rdist = "dM3_nb(mut, g0, g1, rt, J, R, z)",
    discrete = TRUE,
    types = c('value = double(2)',
              'mut = double()',
              'g0  = double()',
              'g1  = double()',
              'rt  = double()',
              'J = integer()',
              'R = integer()',
              'z = double(1)'
              ),
    mixedSizes = FALSE,
    pqAvail = FALSE
  )), verbose = FALSE
)



