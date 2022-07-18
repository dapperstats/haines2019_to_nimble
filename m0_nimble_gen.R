
# multinomial negative binomial mixture model, assuming a vector of observations with indexing based on length vector
#


dM0_nb_vec <- nimbleFunction (
  run = function (x   = double(1),
                  mu  = double(),
                  p   = double(),
                  r   = double(),
                  J_i = integer(1),
                  R   = integer(),
                  log = logical(0, default = 0)) {

  xtot       <- sum(x)
  ptot       <- double(R)
  x_row      <- integer(R)
  x_miss_row <- integer(R)

  for (i in 1:R) {

    spots_in      <- (sum(J_i[1:i]) - J_i[i] + 1):(sum(J_i[1:i]) - J_i[i] + J_i[i])
    ptot[i]       <- 1 - (1 - p)^J_i[i]
    x_row[i]      <- sum(x[spots_in])
    x_miss_row[i] <- x[spots_in] %*% seq(0, J_i[i] - 1)

  }

  x_sumj <- sum(x_miss_row[1:R]) 
  x_logfact     <- sum(lfactorial(x))
  x_row_logfact <- sum(lfactorial(x_row))

  loglik <- (-x_row_logfact)

  for (i in 1:R) {

    spots_in      <- (sum(J_i[1:i]) - J_i[i] + 1):(sum(J_i[1:i]) - J_i[i] + J_i[i])

    x_vec  <- seq(0, J_i[i] - 1)

    term1  <- sum(lgamma(r + x_row[i])) - lgamma(r) 
    term2  <- r * log(r) + x_row[i] * log(mu)
    term3  <- x_row[i] * log(p) + sum(x[spots_in] * x_vec) * log(1 - p)
    term4  <-  -(x_row[i] + r) * log(r + mu * ptot[i])

    loglik <- loglik + term1 + term2 + term3 + term4

  }

  if (log) return(loglik)
  else return(exp(loglik))
  returnType(double())    

})

rM0_nb_vec <- nimbleFunction(
  run = function(n    = integer(),
                 mu   = double(),
                 p    = double(),
                 r    = double(),
                 J_i  = integer(1), 
                 R    = integer()) {

    J_tot <- sum(J_i)    
    

    prob <- double(J_tot + 1 * R)
    retain <- logical(J_tot + 1 * R)

    for (i in 1:R) {

      for (j in 1:J_i[i]) {

        prob[sum(J_i[1:i]) - J_i[i] + i - 1 + j] <- pow(1 - p, j - 1) * p
        retain[sum(J_i[1:i]) - J_i[i] + i - 1 + j] <- TRUE
 
      }
      prob[sum(J_i[1:i]) - J_i[i] + i - 1 + j + 1] <- 1 - sum(prob[(sum(J_i[1:i]) - J_i[i] + i - 1 + 1):(sum(J_i[1:i]) - J_i[i] + i - 1 + J_i[i])])
    }

    ans <- integer(J_tot + 1 * R)
    n <- integer(R)

    for (i in 1:R) {

      n[i] <- rnbinom(n = 1, size = r, mu = mu)

      if (n[i] > 0) {

        ans[(sum(J_i[1:i]) - J_i[i] + i - 1):(sum(J_i[1:i]) - J_i[i] + i - 1 + J_i[i]) + 1] <- rmulti(n = 1, size = n[i], prob = prob[(sum(J_i[1:i]) - J_i[i] + i - 1):(sum(J_i[1:i]) - J_i[i] + i - 1 + J_i[i]) + 1])

      }

    }

    return(ans[retain])
    returnType(double(1))
})


registerDistributions(list(
  dM0_nb_vec = list(
    BUGSdist = "dM0_nb_vec(mu, p, r, J_i, R)",
    Rdist = "dM0_nb_vec(mu, p, r, J_i, R)",
    discrete = TRUE,
    types = c('value = double(1)',
              'mu = double()',
              'p = double()',
              'r = double()',
              'J = integer(1)',
              'R = integer()'
              ),
    mixedSizes = FALSE,
    pqAvail = FALSE
  )), verbose = F
)
