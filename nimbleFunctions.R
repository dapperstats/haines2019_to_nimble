
dM0_nb <- nimbleFunction (
  run = function (x   = integer(2),
                  mu  = double(),
                  p   = double(),
                  r   = double(),
                  J   = integer(),
                  R   = integer(),
                  log = logical(0, default = 0)) {


  ptot <- 1 - (1 - p) ^ J

  x_row      <- integer(R)
  for (i in 1:R) {
    x_row[i] <- sum(x[i, ])
  }

  x_col      <- integer(J)
  for (i in 1:J) {
    x_col[i] <- sum(x[ , i])
  }

  x_tot  <- sum(x)
  x_vec  <- seq(0, J - 1)
  x_sumj <- (x_vec %*% x_col )[1, 1]

  x_logfact   <- sum(lfactorial(x))
  x_r_logfact <- sum(lfactorial(x_row))


  term1   <- sum(lgamma(r + x_row)) - R * lgamma(r) - x_r_logfact
  term2   <- R * r * log(r) + x_tot * log(mu)
  term3   <- x_tot * log(p) + x_sumj * log(1 - p)
  term4   <- -(x_tot + R * r) * log(r + mu * ptot)
  logProb <- term1 + term2 + term3 + term4

  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())  
})



dM0_nb(x   = ymat,
       mu  = log(mut),
       p   = logit(pt),
       r   = log(rt),
       R   = R,
       J   = J,
       log = FALSE)