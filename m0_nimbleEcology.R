####################################################
#                                                  #
# Model 0:                                         #
#          multinomial negative binomial N-mixture #
#                                                  #
#          constant mu                             # 
#          constant p                              #
#                                                  #
#          single site                             #
#                                                  #
####################################################

dNmixture_MNB_s <- nimbleFunction(
    run = function(x   = double(1),
                   mu  = double(),
                   p   = double(),
                   r   = double(),
                   J   = double(),
                   log = integer(0, default = 0)) {

  x_tot <- sum(x)
  x_miss <- sum(x * seq(0, J - 1))

  term1   <- lgamma(r + x_tot) - lgamma(r) - sum(lfactorial(x))
  term2   <- r * log(r) + x_tot * log(mu)
  term3   <- x_tot * log(p) + x_miss * log(1 - p)
  term4   <- -(x_tot + r) * log(r + mu * (1 - (1 - p) ^ J))
  logProb <- term1 + term2 + term3 + term4

  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())  

})


rNmixture_MNB_s <- nimbleFunction(
  run = function(n  = integer(),
                 mu = double(),
                 p  = double(),
                 r  = double(),
                 J  = double()) {


    prob <- numeric(J + 1)
    for (i in 1:(J)) {
      prob[i] <- pow(1 - p, i - 1) * p
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
  dNmixture_MNB_s = list(
    BUGSdist = "dNmixture_MNB_s(mu, p, r, J)",
    Rdist = "dNmixture_MNB_s(mu, p, r, J)",
    discrete = TRUE,
    types = c('value = double(1)',
              'mu = double()',
              'p = double()',
              'r = double()',
              'J = double()'
              ),
    mixedSizes = FALSE,
    pqAvail = FALSE
  )), verbose = FALSE
)