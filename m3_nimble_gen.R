#############################################
# Model 3: removal covariates for capture   #
#############################################

# and 

# not assuming a rectangular survy design but rather 
#   using a vector of observations with indexing based on length vector



dM3_nb_vec <- nimbleFunction (
  run = function (x   = double(1),
                  mut = double(),
                  g0  = double(),
                  g1  = double(),
                  rt  = double(),
                  J_i = integer(),
                  R   = integer(),
                  w   = double(1), 
                  log = logical(0, default = 0)) {

    mu   <- exp(mut)
    pt   <- g0 + g1 * w
    p    <- expit(pt)
    r    <- exp(rt)


  xtot       <- sum(x)
  ptot       <- double(R)
  x_row      <- integer(R)
  x_miss_row <- integer(R)

  for (i in 1:R) {

    spots_in      <- (sum(J_i[1:i]) - J_i[i] + 1):(sum(J_i[1:i]) - J_i[i] + J_i[i])
    x_row[i]      <- sum(x[spots_in])
    x_miss_row[i] <- x[spots_in] %*% seq(0, J_i[i] - 1)

  }

  x_sumj <- sum(x_miss_row[1:R]) 

  x_logfact     <- sum(lfactorial(x))
  
  logProb <- -x_logfact

  for (i in 1:R) {
 
    spots_in      <- (sum(J_i[1:i]) - J_i[i] + 1):(sum(J_i[1:i]) - J_i[i] + J_i[i])

    x_vec  <- seq(0, J_i[i] - 1)
    p_i <- p[spots_in]

    pp_i <- c(0, p_i)
    prob <- double(J_i[i])
    for (j in 1:J_i[i]) {
      prob[j] <- prod(1 - pp_i[1:j]) * p_i[j]
    }
    ptot[i] <- sum(prob)

    term1   <- lgamma(r + x_row[i]) - lgamma(r)
    term2   <- r * log(r) + x_row[i] * log(mu)
    term3   <- sum(x[spots_in] * log(prob))
    term4   <- -(x_row[i] + r) * log(r + mu * ptot[i])
    logProb <- logProb + term1 + term2 + term3 + term4

  }
  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())  
})


rM3_nb_vec <- nimbleFunction(
  run = function(n   = integer(),
                 mut = double(),
                 g0  = double(),
                 g1  = double(),
                 rt  = double(),
                 J_i = integer(), 
                 R   = integer(),
                 w   = double(1)) {

    mu   <- exp(mut)
    pt   <- g0 + g1 * w
    p    <- expit(pt)
    r    <- exp(rt)
    
    J_tot <- sum(J_i)

    prob <- double(J_tot + 1 * R)
    retain <- logical(J_tot + 1 * R)
    ans <- integer(J_tot + 1 * R)


    for (i in 1:R) {

      p_i <- p[(sum(J_i[1:i]) - J_i[i] + 1):sum(J_i[1:i])]

      pp_i <- c(0, p_i)


      for (j in 1:J_i[i]) {

        prob[sum(J_i[1:i]) - J_i[i] + i - 1 + j] <- prod(1 - pp_i[1:j]) * p_i[j]
        retain[sum(J_i[1:i]) - J_i[i] + i - 1 + j] <- TRUE
 
      }
      prob[sum(J_i[1:i]) - J_i[i] + i - 1 + j + 1] <- 1 - sum(prob[(sum(J_i[1:i]) - J_i[i] + i - 1 + 1):(sum(J_i[1:i]) - J_i[i] + i - 1 + J_i[i])])
    }

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
  dM3_nb_vec = list(
    BUGSdist = "dM3_nb_vec(mut, g0, g1, rt, J_i, R, z)",
    Rdist = "dM3_nb_vec(mut, g0, g1, rt, J_i, R, z)",
    discrete = TRUE,
    types = c('value = double(1)',
              'mut = double()',
              'g0  = double()',
              'g1  = double()',
              'rt  = double()',
              'J_i = double()',
              'R = integer()',
              'z = double(1)'
              ),
    mixedSizes = FALSE,
    pqAvail = FALSE
  )), verbose = FALSE
)



