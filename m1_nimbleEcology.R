####################################################
#                                                  #
# Model 1:                                         #
#          multinomial negative binomial N-mixture #
#                                                  #
#          covariates for mu                       # 
#          constant p                              #
#                                                  #
#          single site                             #
#                                                  #
####################################################

dNmixture_MNB_sitecovar_s <- nimbleFunction(
    run = function(x   = double(1),
                   b0  = double(),
                   b1  = double(),
                   p   = double(),
                   r   = double(),
                   J   = integer(),
                   z   = double(),
                   log = integer(0, default = 0)) {

  zb  <- log(b0) + b1 * z # "log mu"
  mu  <- exp(zb)

  x_tot <- sum(x)
  x_miss <- sum(x * seq(0, J - 1))

  term1   <- lgamma(r + x_tot) - lgamma(r) - sum(lfactorial(x))
  term2   <- r * log(r) + x_tot * zb
  term3   <- x_tot * log(p) + x_miss * log(1 - p)
  term4   <- -(x_tot + r) * log(r + mu * (1 - (1 - p) ^ J))
  logProb <- term1 + term2 + term3 + term4

  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())  

})


rNmixture_MNB_sitecovar_s <- nimbleFunction(
  run = function(n  = integer(),
                 b0  = double(),
                 b1  = double(),
                 p  = double(),
                 r  = double(),
                 J  = integer(),
                 z   = double()) {


    zb  <- log(b0) + b1 * z
    mu  <- exp(zb)
        
    prob <- double(J + 1)
    for (i in 1:(J)) {
      prob[i] <- pow(1 - p, i - 1) * p
    }
    prob[J + 1] <- 1 - sum(prob[1:J])

    ans <- integer(J + 1)
    n <- rnbinom(n = 1, size = r, mu = mu)
    if (n > 0) {
      ans <- rmulti(n = 1, size = n, prob = prob)
    }

    return(ans[1:J])
    returnType(integer(1))
})

registerDistributions(list(
  dNmixture_MNB_sitecovar_s = list(
    BUGSdist = "dNmixture_MNB_sitecovar_s(b1, b0, p, r, J, z)",
    Rdist = "dNmixture_MNB_sitecovar_s(b1, b0, p, r, J, z)",
    discrete = TRUE,
    types = c('value = double(1)',
              'b0 = double()',
              'b1 = double()',
              'p = double()',
              'r = double()',
              'J = integer()',
              'z = double()'
              ),
    mixedSizes = FALSE,
    pqAvail = FALSE
  )), verbose = FALSE
)