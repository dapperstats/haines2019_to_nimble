####################################################
#                                                  #
# Model 3:                                         #
#          multinomial negative binomial N-mixture #
#                                                  #
#          constant mu                             # 
#          dynamic  p                              #
#                                                  #
#          single site                             #
#                                                  #
####################################################

dNmixture_MNB_v <- nimbleFunction(
    run = function(x   = double(1),
                   mu  = double(),
                   p   = double(1),
                   r   = double(),
                   J   = double(),
                   log = integer(0, default = 0)) {

  x_tot <- sum(x)
  x_miss <- sum(x * seq(0, J - 1))
  pp <- c(0, p)
    prob <- numeric(J)
    for (j in 1:J) {
      prob[j] <- prod(1 - pp[1:j]) * p[j]
    }
    ptot <- sum(prob)

  term1   <- lgamma(r + x_tot) - lgamma(r) - sum(lfactorial(x))
  term2   <- r * log(r) + x_tot * log(mu)
  term3   <- sum(x * log(prob))
  term4   <- -(x_tot + r) * log(r + mu * ptot)
  logProb <- term1 + term2 + term3 + term4

  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())  

})


rNmixture_MNB_v <- nimbleFunction(
  run = function(n  = integer(),
                 mu = double(),
                 p  = double(1),
                 r  = double(),
                 J  = double()) {


    pp <- c(0, p)
    prob <- numeric(J + 1)
    for (j in 1:J) {
      prob[j] <- prod(1 - pp[1:j]) * p[j]
    }
    prob[J + 1] <- 1 - sum(prob[1:J])

    ans <- numeric(J + 1)
    n <- rnbinom(n = 1, size = r, prob = 1/(1 + (1/r) * mu))
    if (n > 0) {
      ans <- rmulti(n = 1, size = n, prob = prob)
    }

    return(ans[1:J])
    returnType(double(1))
})

registerDistributions(list(
  dNmixture_MNB_v = list(
    BUGSdist = "dNmixture_MNB_v(mu, p, r, J)",
    Rdist = "dNmixture_MNB_v(mu, p, r, J)",
    discrete = TRUE,
    types = c('value = double(1)',
              'mu = double()',
              'p = double(1)',
              'r = double()',
              'J = double()'
              ),
    mixedSizes = FALSE,
    pqAvail = FALSE
  )), verbose = FALSE
)