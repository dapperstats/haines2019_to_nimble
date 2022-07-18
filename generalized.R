
# multinomial negative binomial mixture model, assuming a vector of observations with indexing based on length vector
#

# working here! the randomizer function works great
# need to isolate the "colSum" to figure out how to update the code

dM0_nb <- nimbleFunction (
  run = function (x   = double(1),
                  mu  = double(),
                  p   = double(),
                  r   = double(),
                  J_i = integer(1),
                  R   = integer(),
                  log = logical(0, default = 0)) {


  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())    

})

rM0_nb_vec <- nimbleFunction(
  run = function(n    = integer(),
                 mu   = double(),
                 p    = double(),
                 r    = double(),
                 J_i  = integer(1), 
                 R    = integer()) {

    J_tot <- sum(J_i)    
    

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

      n[i] <- rnbinom(n = 1, size = r, mu = mu)

      if (n[i] > 0) {

        ans[(sum(J_i[1:i]) - J_i[i] + i - 1):(sum(J_i[1:i]) - J_i[i] + i - 1 + J_i[i]) + 1] <- rmulti(n = 1, size = n[i], prob = prob[(sum(J_i[1:i]) - J_i[i] + i - 1):(sum(J_i[1:i]) - J_i[i] + i - 1 + J_i[i]) + 1])

      }

    }

    return(ans[retain])
    returnType(double(1))
})
