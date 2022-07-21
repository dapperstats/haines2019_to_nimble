#
# converting the basic m0 from haines into a carcass model
#


# this model still assumed rectangular data 

# alpha = arrival rate
#  alphat = transformed (log-scale) arrival rate
# 
# omega = observation probability
#  omegat = transformed (logit-scale) observation probability

# D = day of J (samples)

#   WORKING IN HERE BUT THIS IS NOT QUITE RIGHT OR DONE YET!

dM0_nb <- nimbleFunction (
  run = function (x     = double(2),
                 alphat = double(),
                 omegat = double(),
                 rt     = double(),
                 J      = integer(), 
                 R      = integer(), 
                 D      = integer(1),
                 log    = logical(0, default = 0)) {

    alpha <- exp(alphat)
    omega <- expit(omegat)
    r     <- exp(rt)
    
    DD <- c(0, D)

    alphas <- double(J)
    for (j in 1:J) {
      alphas[j] <- alpha * (DD[j+1] - DD[j])
    }

   omegas <- rep(omega, J)
   pp <- c(0, omegas)


    prob <- double(J)
    for (j in 1:J) {
      prob[j] <- prod(1 - pp[1:j]) * omegas[j]
    }
    ptot <- sum(prob)

  x_row      <- integer(R)
  x_miss_row <- integer(R)

  for (i in 1:R) {
    x_row[i]      <- sum(x[i, ])
    x_miss_row[i] <- x[i, ] %*% seq(0, J - 1)
  }

  x_col <- integer(J)
  for (j in 1:J) {
    x_col[j] <- sum(x[,j])
  }

  x_tot  <- sum(x)
  x_vec  <- seq(0, J - 1)
  x_sumj <- sum(x_miss_row[1:R]) 

  x_logfact     <- sum(lfactorial(x))
  x_row_logfact <- sum(lfactorial(x_row))

  term1   <- sum(lgamma(r + x_row)) - R * lgamma(r) - x_logfact
  term2   <- R * r * log(r) +  x_col %*% log(alphas) 
  term3   <- x_tot * log(omega) + x_sumj * log(1 - omega)
  term4   <- -(x_tot + R * r) * log(r + alphas %*% prob)
  logProb <- term1 + term2 + term3 + term4

  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())  
})




rM0_nb_carcass <- nimbleFunction(
  run = function(n      = integer(),
                 alphat = double(),
                 omegat = double(),
                 rt     = double(),
                 J      = integer(), 
                 R      = integer(), 
                 D      = integer(1)) {

    alpha <- exp(alphat)
    omega <- expit(omegat)
    r     <- exp(rt)
    
    DD <- c(0, D)

    prob   <- double(J + 1)
    alphas <- double(J)
    for (j in 1:J) {
      prob[j]   <- pow(1 - omega, j - 1) * omega
      alphas[j] <- alpha * (DD[j+1] - DD[j])
    }
    prob[J + 1] <- 1 - sum(prob[1:J])

    alpha_tot <- sum(alphas)

    ans <- matrix(0, nrow = R, ncol = J + 1)
    for (i in 1:R) {
      n <- rnbinom(n = 1, size = r, mu = alpha_tot)
      if (n > 0) {
        ans[i, ] <- rmulti(n = 1, size = n, prob = prob)
      }
    }

    return(ans[ , 1:J])
    returnType(integer(2))
})
