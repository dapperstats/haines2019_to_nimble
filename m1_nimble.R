
# multinomial negative binomial mixture model
#   assuming a rectangular  matrix of observations
#   but generalizing to allow more flexibles

# this formulation assumes that the whole matrix of observations is dependent

# i'm not certain the script in the ms is correct for this model
#   it has an extra function that isnt used and there is a missing component in the calculations
#       the logfact is not in the initial value of logProb or in term1
# i've added the term in here 


dM1_nb <- nimbleFunction (
  run = function (x   = double(2),
                  b0  = double(),
                  b1  = double(),
                  pt  = double(),
                  rt  = double(),
                  J   = integer(),
                  R   = integer(),
                  z   = double(1), 
                  log = logical(0, default = 0)) {

  mut <- b0 + b1 * z
  mu   <- exp(mut)
  p    <- expit(pt)
  r    <- exp(rt)

  ptot <- 1 - (1 - p)^J

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

    term1   <- lgamma(r + x_row[i]) - lgamma(r)
    term2   <- r * log(r) + x_row[i] * mut[i]
    term3   <- x_row[i] * log(p) + x_miss_row[i] * log(1 - p)
    term4   <- -(x_row[i] + r) * log(r + mu[i] * ptot)
    logProb <- logProb + term1 + term2 + term3 + term4

  }
  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())  
})


rM1_nb <- nimbleFunction(
  run = function(n  = integer(),
                 b0 = double(),
                 b1 = double(),
                 pt = double(),
                 rt = double(),
                 J  = integer(), 
                 R  = integer(),
                 z  = double(1)) {

    mut <- b0 + b1 * z
    mu   <- exp(mut)
    p    <- expit(pt)
    r    <- exp(rt)
    
    prob <- double(J + 1)
    for (i in 1:(J)) {
      prob[i] <- pow(1 - p, i - 1) * p
    }
    prob[J + 1] <- 1 - sum(prob[1:J])

    ans <- matrix(0, nrow = R, ncol = J + 1)

    for (i in 1:R) {
      n <- rnbinom(n = 1, size = r, mu = mu[i])
      if (n > 0) {
        ans[i, ] <- rmulti(n = 1, size = n, prob = prob)
      }
    }

    return(ans[ , 1:J])
    returnType(integer(2))
})

registerDistributions(list(
  dM1_nb = list(
    BUGSdist = "dM1_nb(b0, b1, pt, rt, J, R, z)",
    Rdist = "dM1_nb(b0, b1, pt, rt, J, R, z)",
    discrete = TRUE,
    types = c('value = double(2)',
              'b0  = double()',
              'b1  = double()',
              'pt  = double()',
              'rt  = double()',
              'J = integer()',
              'R = integer()',
              'z = double(1)'
              ),
    mixedSizes = FALSE,
    pqAvail = FALSE
  )), verbose = FALSE
)



