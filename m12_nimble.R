# combining

##########################################
# Model 1: site covariates for abundance #
##########################################

# and

##########################################
# Model 2: site covariates for capture   #
##########################################


dM12_nb <- nimbleFunction (
  run = function (x   = double(2),
                  b0  = double(),
                  b1  = double(),
                  g0  = double(),
                  g1  = double(),
                  rt  = double(),
                  J   = integer(),
                  R   = integer(),
                  z   = double(1),
                  w   = double(1), 
                  log = logical(0, default = 0)) {


  mu   <- exp(b0 + b1 * z)
  pt   <- g0 + g1 * w
  p    <- expit(pt)

  r    <- exp(rt)

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

    prob <- numeric(J)
    for (j in 1:J) {
      prob[j] <- pow(1 - p[i], j - 1) * p[i]
    }
    ptot <- sum(prob)

    term1   <- lgamma(r + x_row[i]) - lgamma(r)
    term2   <- r * log(r) + x_row[i] * log(mu[i])
    term3   <- sum(ymat[i,] * log(prob))
    term4   <- -(x_row[i] + r) * log(r + mu[i] * ptot)
    logProb <- logProb + term1 + term2 + term3 + term4

  }
  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())  
})


rM12_nb <- nimbleFunction(
  run = function(n  = integer(),
                 b0 = double(),
                 b1 = double(),
                 g0  = double(),
                 g1  = double(),
                 rt = double(),
                 J  = integer(), 
                 R  = integer(),
                 z  = double(1),
                 w  = double(1)) {

    mut <- b0 + b1 * z
    mu   <- exp(mut)
    pt   <- g0 + g1 * w
    p    <- expit(pt)

    r    <- exp(rt)
    
    ans <- matrix(0, nrow = R, ncol = J + 1)

    for (i in 1:R) {
      n <- rnbinom(n = 1, size = r, mu = mu[i])
      if (n > 0) {

        prob <- double(J + 1)
        for (i in 1:J) {
          prob[i] <- pow(1 - p[i], j - 1) * p[i]
        }
        prob[J + 1] <- 1 - sum(prob[1:J])

        ans[i, ] <- rmulti(n = 1, size = n, prob = prob)
      }
    }

    return(ans[ , 1:J])
    returnType(integer(2))
})

registerDistributions(list(
  dM12_nb = list(
    BUGSdist = "dM12_nb(b0, b1, g0, g1, rt, J, R, z, w)",
    Rdist = "dM12_nb(b0, b1, g0, g1, rt, J, R, z, w)",
    discrete = TRUE,
    types = c('value = double(2)',
              'b0  = double()',
              'b1  = double()',
              'g0  = double()',
              'g1  = double()',
              'rt  = double()',
              'J = integer()',
              'R = integer()',
              'z = double(1)',
              'w = double(1)'
              ),
    mixedSizes = FALSE,
    pqAvail = FALSE
  )), verbose = FALSE
)



