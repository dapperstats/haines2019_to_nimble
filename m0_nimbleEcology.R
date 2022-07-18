####################################################
# Model 0:                                         #
#          multinomial negative binomial N-mixture #
#          constant mu                             # 
#           constant p                             #
####################################################

dNmixture_MNB_s <- nimbleFunction(
    run = function(x      = double(1),
                   mu  = double(),
                   p   = double(),
                   r   = double(),
                   J   = integer(),
                   log = integer(0, default = 0)) {


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