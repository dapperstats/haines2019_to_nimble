##########################################
# Model 1: site covariates for abundance #
##########################################

# and 

# not assuming a rectangular survy design but rather 
#   using a vector of observations with indexing based on length vector


dM1_nb_vec <- nimbleFunction (
  run = function (x   = double(1),
                  b0  = double(),
                  b1  = double(),
                  pt  = double(),
                  rt  = double(),
                  J_i = integer(1),
                  R   = integer(),
                  z   = double(1), 
                  log = logical(0, default = 0)) {

  mut <- b0 + b1 * z
  mu   <- exp(mut)
  p    <- expit(pt)
  r    <- exp(rt)

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

  loglik <- (-x_logfact) 

  for (i in 1:R) {

    spots_in      <- (sum(J_i[1:i]) - J_i[i] + 1):(sum(J_i[1:i]) - J_i[i] + J_i[i])

    x_vec  <- seq(0, J_i[i] - 1)

    term1  <- sum(lgamma(r + x_row[i])) - lgamma(r) 
    term2  <- r * log(r) + x_row[i] * log(mu[i])
    term3  <- x_row[i] * log(p) + sum(x[spots_in] * x_vec) * log(1 - p)
    term4  <-  -(x_row[i] + r) * log(r + mu[i] * ptot[i])

    loglik <- loglik + term1 + term2 + term3 + term4

  }

  if (log) return(loglik)
  else return(exp(loglik))
  returnType(double())      
})


rM1_nb_vec <- nimbleFunction(
  run = function(n   = integer(),
                 b0  = double(),
                 b1  = double(),
                 pt  = double(),
                 rt  = double(),
                 J_i = double(), 
                 R   = integer(),
                 z   = double(1)) {

    mut <- b0 + b1 * z
    mu   <- exp(mut)
    p    <- expit(pt)
    r    <- exp(rt)
    
    J_tot <- sum(J_i)

    prob <- double(J_tot + 1 * R)

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

      n[i] <- rnbinom(n = 1, size = r, mu = mu[i])

      if (n[i] > 0) {

        ans[(sum(J_i[1:i]) - J_i[i] + i - 1):(sum(J_i[1:i]) - J_i[i] + i - 1 + J_i[i]) + 1] <- rmulti(n = 1, size = n[i], prob = prob[(sum(J_i[1:i]) - J_i[i] + i - 1):(sum(J_i[1:i]) - J_i[i] + i - 1 + J_i[i]) + 1])

      }

    }

    return(ans[retain])
    returnType(double(1))
})

registerDistributions(list(
  dM1_nb_vec = list(
    BUGSdist = "dM1_nb_vec(b0, b1, pt, rt, J_i, R, z)",
    Rdist = "dM1_nb_vec(b0, b1, pt, rt, J_i, R, z)",
    discrete = TRUE,
    types = c('value = double(1)',
              'b0  = double()',
              'b1  = double()',
              'pt  = double()',
              'rt  = double()',
              'J_i = double(1)',
              'R = integer()',
              'z = double(1)'
              ),
    mixedSizes = FALSE,
    pqAvail = FALSE
  )), verbose = FALSE
)



