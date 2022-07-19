
# multinomial negative binomial mixture model
#   assuming a rectangular  matrix of observations
#   but generalizing to allow more flexibles

# this formulation assumes that the whole matrix of observations is dependent
# the term3 calculation is represented differently here, as in the ms
# this is the identical calculation represented like in the other models
# term3 <-  x_row[i] * log(p[i]) + sum(x[i, ] * x_vec) * log(1 - p[i]) 
# it's just harder to see the variation over j

dM2_nb <- nimbleFunction (
  run = function (x   = double(2),
                  mut = double(),
                  g0  = double(),
                  g1  = double(),
                  rt  = double(),
                  J   = integer(),
                  R   = integer(),
                  w   = double(1), 
                  log = logical(0, default = 0)) {

    mu   <- exp(mut)
    pt   <- g0 + g1 * w
    p    <- expit(pt)
    r    <- exp(rt)

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
  
  logProb <- -x_logfact

  for (i in 1:R) {
 
    ptot <- 1 - (1 - p[i])^J
   
    prob <- double(J)
    for (j in 1:J) {
      prob[j] <- pow(1 - p[i], j - 1) * p[i]
    }

    term1   <- lgamma(r + x_row[i]) - lgamma(r)
    term2   <- r * log(r) + x_row[i] * log(mu)
    term3   <- sum(x[i, ] * log(prob))
    term4   <- -(x_row[i] + r) * log(r + mu * ptot)
    logProb <- logProb + term1 + term2 + term3 + term4

  }
  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())  
})


rM2_nb <- nimbleFunction(
  run = function(n   = integer(),
                 mut = double(),
                 g0  = double(),
                 g1  = double(),
                 rt  = double(),
                 J   = integer(), 
                 R   = integer(),
                 w   = double(1)) {

    mu   <- exp(mut)
    pt   <- g0 + g1 * w
    p    <- expit(pt)
    r    <- exp(rt)
    
    ans <- matrix(0, nrow = R, ncol = J + 1)

    for (i in 1:R) {
      n <- rnbinom(n = 1, size = r, mu = mu)
      if (n > 0) {
        prob <- double(J + 1)
        for (j in 1:J) {
          prob[j] <- pow(1 - p[i], j - 1) * p[i]
        }
        prob[J + 1] <- 1 - sum(prob[1:J])


        ans[i, ] <- rmulti(n = 1, size = n, prob = prob)
      }
    }

    return(ans[ , 1:J])
    returnType(integer(2))
})

registerDistributions(list(
  dM2_nb = list(
    BUGSdist = "dM2_nb(mut, g0, g1, rt, J, R, z)",
    Rdist = "dM2_nb(mut, g0, g1, rt, J, R, z)",
    discrete = TRUE,
    types = c('value = double(2)',
              'mut = double()',
              'g0  = double()',
              'g1  = double()',
              'rt  = double()',
              'J = integer()',
              'R = integer()',
              'z = double(1)'
              ),
    mixedSizes = FALSE,
    pqAvail = FALSE
  )), verbose = FALSE
)



