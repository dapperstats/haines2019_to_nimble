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
