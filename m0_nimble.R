
# multinomial negative binomial mixture model
#   assuming a rectangular  matrix of observations
#   but generalizing to allow more flexibles

# this formulation assumes that the whole matrix of observations is dependent

dM0_nb <- nimbleFunction (
  run = function (x   = double(2),
                  mu  = double(),
                  p   = double(),
                  r   = double(),
                  J   = integer(),
                  R   = integer(),
                  log = logical(0, default = 0)) {


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
  x_row_logfact <- sum(lfactorial(x_row))

  term1   <- sum(lgamma(r + x_row)) - R * lgamma(r) - x_row_logfact
  term2   <- R * r * log(r) + x_tot * log(mu)
  term3   <- x_tot * log(p) + x_sumj * log(1 - p)
  term4   <- -(x_tot + R * r) * log(r + mu * ptot)
  logProb <- term1 + term2 + term3 + term4

  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())  
})


rM0_nb <- nimbleFunction(
  run = function(n  = integer(),
                 mu = double(),
                 p  = double(),
                 r  = double(),
                 J  = integer(), 
                 R  = integer()) {

    
    prob <- double(J + 1)
    for (i in 1:(J)) {
      prob[i] <- pow(1 - p, i - 1) * p
    }
    prob[J + 1] <- 1 - sum(prob[1:J])

    ans <- matrix(0, nrow = R, ncol = J + 1)
    for (i in 1:R) {
      n <- rnbinom(n = 1, size = r, mu = mu)
      if (n > 0) {
        ans[i, ] <- rmulti(n = 1, size = n, prob = prob)
      }
    }

    return(ans[ , 1:J])
    returnType(integer(2))
})

registerDistributions(list(
  dM0_nb = list(
    BUGSdist = "dM0_nb(mu, p, r, J, R)",
    Rdist = "dM0_nb(mu, p, r, J, R)",
    discrete = TRUE,
    types = c('value = double(2)',
              'mu = double()',
              'p = double()',
              'r = double()',
              'J = integer()',
              'R = integer()'
              ),
    mixedSizes = FALSE,
    pqAvail = FALSE
  )), verbose = F
)

