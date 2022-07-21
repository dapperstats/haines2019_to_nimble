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

# and 

# not assuming a rectangular survy design but rather 
#   using a vector of observations with indexing based on length vector

dM123_nb_vec <- nimbleFunction (
  run = function (x   = double(1),
                  b0  = double(),
                  b1  = double(),
                  g0  = double(),
                  g1  = double(),
                  g2  = double(),
                  rt  = double(),
                  J_i = double(1),
                  R   = integer(),
                  z   = double(1),
                  wR  = double(1), 
                  wJ  = double(1), 
                  log = logical(0, default = 0)) {


  mu   <- exp(b0 + b1 * z)


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

    p    <- expit(g0 + g1 * wR[i] + g2 * wJ)
    pp <- c(0, p)

    prob <- numeric(J)
    for (j in 1:J) {
      prob[j] <- prod(1 - pp[1:j])*p[j]
    }
    ptot <- sum(prob)

    term1   <- lgamma(r + x_row[i]) - lgamma(r)
    term2   <- r * log(r) + x_row[i] * log(mu[i])
    term3   <- sum(x[i,] * log(prob))
    term4   <- -(x_row[i] + r) * log(r + mu[i] * ptot)
    logProb <- logProb + term1 + term2 + term3 + term4

  }
  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())  
})


rM123_nb_vec <- nimbleFunction(
  run = function(n  = integer(),
                 b0 = double(),
                 b1 = double(),
                 g0  = double(),
                 g1  = double(),
                 g2  = double(),
                 rt  = double(),
                 J_i = integer(), 
                 R   = integer(),
                 z  = double(1),
                 wR = double(1),
                 wJ = double(1)) {


    mu   <- exp(b0 + b1 * z)

    r    <- exp(rt)
    
    ans <- matrix(0, nrow = R, ncol = J + 1)

    for (i in 1:R) {
      n <- rnbinom(n = 1, size = r, mu = mu[i])
      if (n > 0) {

        p    <- expit(g0 + g1 * wR[i] + g2 * wJ)
        pp <- c(0, p)

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
  dM123_nb_vec = list(
    BUGSdist = "dM123_nb_vec(b0, b1, g0, g1, g2, rt, J_i, R, z, wR, wJ)",
    Rdist = "dM123_nb_vec(b0, b1, g0, g1, g2, rt, J_i, R, z, wR, wJ)",
    discrete = TRUE,
    types = c('value = double(1)',
              'b0  = double()',
              'b1  = double()',
              'g0  = double()',
              'g1  = double()',
              'g2  = double()',
              'rt  = double()',
              'J_i = double(1)',
              'R = integer()',
              'z = double(1)',
              'wR = double(1)',
              'wJ = double(1)'
              ),
    mixedSizes = FALSE,
    pqAvail = FALSE
  )), verbose = FALSE
)



