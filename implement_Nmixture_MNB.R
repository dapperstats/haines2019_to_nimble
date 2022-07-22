#
#  This script details the use of two Multinomial - Negative Binomial mixture models in nimble
#
#   the distributions follow the approach of function naming for nimbleEcology
#
#   Nmixture_MNB_s  is a multinomial with scalar p 
#   Nmixture_MNB_v  is a multinomial with vector p
#
#  functions are parameterized in terms of mean (mu), scale (r) and visit-specific detection probability (p)
#  there are J visits in the multinomial
#
#  the distributions are registered in their script and the models (and their samplers) are compiled below
#
rm(list = ls())

library(nimble)

logit <- function (x) {
   log(x / (1 - x))
}
expit <- function (x) { 
  exp(x) / (1 + exp(x))
}

source("Nmixture_MNB.R")

set.seed(1312)

mu <- 28.05
p <- 0.61
r <- 11.4
rt <- log(r)
b0 <- 3.077
b1 <- 0.101
g0 <- 0.5
g1 <- 0.125
g2 <- 0.125

R   <- 20
J_i <- rep(3:5, length.out = R)
J_tot <- sum(J_i)

spot1 <- integer(R)
spot2 <- integer(R)
for (i in 1:R) {
  spot1[i] <- (sum(J_i[1:i]) - J_i[i] + 1)
  spot2[i] <- (sum(J_i[1:i]))
}


# covariates

z  <- rnorm(R, 0, 1)
wR_draw <- rnorm(R, 0, 1) 
wR      <- rep(wR_draw, J_i)
wJ_draw <- rnorm(max(J_i), 0, 1)
wJ      <- wJ_draw[sequence(J_i)]

mu_i  <- exp(b0 + b1 * z)
p_s_i <- expit(g0 + g1 * wR_draw)
p_v_j <- expit(g0 + g1 * wR + g2 * wJ)

# generate from scalar p
ys <- NULL
for (i in 1:R) {
  ys <- c(ys, rNmixture_MNB_s(n = 1, mu = mu_i[i], p = p_s_i[i], r = r, J = J_i[i]))
}

# generate from vector p
yv <- NULL
for (i in 1:R) {
  spots_in <- (sum(J_i[1:i]) - J_i[i] + 1):sum(J_i[1:i])
  yv <- c(yv, rNmixture_MNB_v(n = 1, mu = mu_i[i], p = p_v_j[spots_in], r = r, J = J_i[i]))
}


nc_s <- nimbleCode({

  rt ~ dnorm(0, 0.1)

  b0 ~ dnorm(0, 0.5)
  b1 ~ dnorm(0, 0.1)

  g0 ~ dnorm(0, 0.5)
  g1 ~ dnorm(0, 0.1)


  r <- exp(rt)
  
  for (i in 1:R) {
    mut_i[i] <- b0 + b1 * z[i] 
    mu_i[i]  <- exp(mut_i[i])
    pt_i[i] <- g0 + g1 * wR[i]
    p_i[i] <- expit(pt_i[i])

    x[spot1[i]:spot2[i]] ~ dNmixture_MNB_s(mu = mu_i[i], p = p_i[i], r = r, J = J_i[i])
  }
})


nc_v <- nimbleCode({

  rt ~ dnorm(0, 0.1)

  b0 ~ dnorm(0, 0.5)
  b1 ~ dnorm(0, 0.1)

  g0 ~ dnorm(0, 0.5)
  g1 ~ dnorm(0, 0.1)
  g2 ~ dnorm(0, 0.1)

  r <- exp(rt)
  
  pt_j[1:J_tot] <- g0 + g1 * wR[1:J_tot] + g2 * wJ[1:J_tot]
  p_j[1:J_tot] <- expit(pt_j[1:J_tot])

  for (i in 1:R) {
    mut_i[i] <- b0 + b1 * z[i] 
    mu_i[i]  <- exp(mut_i[i])
    x[spot1[i]:spot2[i]] ~ dNmixture_MNB_v(mu = mu_i[i], p = p_j[spot1[i]:spot2[i]], r = r, J = J_i[i])
  }
})








# scalar data, scalar model

nmix_ss <- nimbleModel(nc_s,
                      constants = list(J_i = J_i, R = R, spot1 = spot1, spot2 = spot2, J_tot = J_tot),
                      data = list(x = ys, z = z, wR = wR),
                      inits = list(b0 = b0,
                                   b1 = b1,
                                   g0 = g0,
                                   g1 = g1,
                                   rt = rt))
nmix_ss$calculate()



# scalar data, vector model

nmix_sv <- nimbleModel(nc_v,
                      constants = list(J_i = J_i, R = R, spot1 = spot1, spot2 = spot2, J_tot = J_tot),
                      data = list(x = ys, z = z, wR = wR, wJ = wJ),
                      inits = list(b0 = b0,
                                   b1 = b1,
                                   g0 = g0,
                                   g1 = g1,
                                   g2 = g2,
                                   rt = rt))
nmix_sv$calculate()





# vector data, scalar model

nmix_vs <- nimbleModel(nc_s,
                      constants = list(J_i = J_i, R = R, spot1 = spot1, spot2 = spot2, J_tot = J_tot),
                      data = list(x = yv, z = z, wR = wR),
                      inits = list(b0 = b0,
                                   b1 = b1,
                                   g0 = g0,
                                   g1 = g1,
                                   rt = rt))
nmix_vs$calculate()



# vector data, vector model

nmix_vv <- nimbleModel(nc_v,
                      constants = list(J_i = J_i, R = R, spot1 = spot1, spot2 = spot2, J_tot = J_tot),
                      data = list(x = yv, z = z, wR = wR, wJ = wJ),
                      inits = list(b0 = b0,
                                   b1 = b1,
                                   g0 = g0,
                                   g1 = g1,
                                   g2 = g2,
                                   rt = rt))
nmix_vv$calculate()



nmix_ss$calculate()
nmix_sv$calculate()
nmix_vs$calculate()
nmix_vv$calculate()


# compile and sample the models 10,000 times


compileNimble(nmix_ss)
nmix_ssMCMC <- buildMCMC(nmix_ss)
cnmix_ssMCMC <- compileNimble(nmix_ssMCMC, project = nmix_ss)
samples_ss <- runMCMC(cnmix_ssMCMC, niter = 10000)

samples_ss
plot(samples_ss[, "b0"], type = "l")

compileNimble(nmix_vv)
nmix_vvMCMC <- buildMCMC(nmix_vv)
cnmix_vvMCMC <- compileNimble(nmix_vvMCMC, project = nmix_vv)
samples_vv <- runMCMC(cnmix_vvMCMC, niter = 10000)

samples_vv
plot(samples_vv[, "b0"], type = "l")



