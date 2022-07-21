rm(list = ls())

library(nimble)

logit <- function (x) {
   log(x / (1 - x))
}
expit <- function (x) { 
  exp(x) / (1 + exp(x))
}

###################################
#           M0                    #
###################################

source("m0_haines.R")
source("m0_nimble.R")
source("m0_nimble_gen.R")
source("m0_nimbleEcology.R")

R  <- 20
J  <- 3
J_i <- rep(J, R)
set.seed(125)

p   <- 0.61
mu  <- 28.048
r   <- 11.415

pt  <- logit(p)
mut <- log(mu)
rt  <- log(r)


param <- c(mut, pt, rt)

ymat  <- gM0nbgen(param, J, R)
ymatv <- as.numeric(t(ymat))

gM0nb(param, ymat, J, R)
dM0_nb(ymat, mut, pt, rt, J, R, TRUE)
dM0_nb_vec(ymatv, mut, pt, rt, J_i, R, TRUE)


dNmixture_MNB_s(ymat[1, ], mut, pt, rt, J, TRUE)
dM0_nb(ymat[1, , drop = FALSE], mut, pt, rt, J, 1, TRUE)
gM0nb(param, ymat[1, , drop = FALSE], J, 1)
dM0_nb_vec(ymat[1, ], mut, pt, rt, J, 1, TRUE)


set.seed(123)
rM0_nb(1, mut, pt, rt, J, R)
set.seed(123)
rNmixture_MNB_s(1, mut, pt, rt, J)
set.seed(123)
rM0_nb_vec(1, mut, pt, rt, J_i, R)
set.seed(123)
gM0nbgen(param, J, R)


J_i_uneven <- rep(2:4, length.out = R)
ymatv_uneven <- rM0_nb_vec(1, mut, pt, rt, J_i_uneven, R)
dM0_nb_vec(ymatv_uneven, mut, pt, rt, J_i_uneven, R, TRUE)



nc <- nimbleCode({
   x[1:J] ~ dNmixture_MNB_s(mut = mut, pt = logit(p), rt = rt, J = J)

 })


nmix <- nimbleModel(nc,
                    constants = list(J = J),
                    data = list(x = ymatv[1:J]),
                    inits = list(mut = mut,
                                 p = p,
                                 rt = rt))
nmix$calculate()

dNmixture_MNB_s(ymatv[1:J], mut, pt, rt, J, log = 1)
dM0_nb_vec(ymat[1, ], mut, pt, rt, J, 1, TRUE)



nc <- nimbleCode({


   for (i in 1:R) {

     x[spot1[i]:spot2[i]] ~ dNmixture_MNB_s(mut = log(mu), pt = logit(p), rt = rt, J = J_i[i])
   }
 })

spot1 <- integer(R)
spot2 <- integer(R)
for (i in 1:R) {
  spot1[i] <- (sum(J_i[1:i]) - J_i[i] + 1)
  spot2[i] <- (sum(J_i[1:i]))
}

J_tot <- sum(J_i)

nmix <- nimbleModel(nc,
                    constants = list(J_i = J_i, R = R, spot1 = spot1, spot2 = spot2),
                    data = list(x = ymatv),
                    inits = list(mu = mu,
                                 p = p,
                                 rt = rt))

nmix$calculate()
gM0nb(param, ymat, J, R)
dM0_nb(ymat, mut, pt, rt, J, R, TRUE)
dM0_nb_vec(ymatv, mut, pt, rt, J_i, R, TRUE)

     spot1_uneven <- integer(R)
     spot2_uneven <- integer(R)
for (i in 1:R) {
     spot1_uneven[i] <- (sum(J_i_uneven[1:i]) - J_i_uneven[i] + 1)
     spot2_uneven[i] <- (sum(J_i_uneven[1:i]))
}
J_tot_uneven <- sum(J_i_uneven)

nmix <- nimbleModel(nc,
                    constants = list(J_i = J_i_uneven, R = R, spot1 = spot1_uneven, spot2 = spot2_uneven),
                    data = list(x = ymatv_uneven),
                    inits = list(mu = mu,
                                 p = p,
                                 rt = rt))

nmix$calculate()
dM0_nb_vec(ymatv_uneven, mut, pt, rt, J_i_uneven, R, TRUE)



###################################
#           M1                    #
###################################

source("m1_haines.R")
source("m1_nimble.R")
source("m1_nimble_gen.R")


b0 <- 3.077
b1 <- 0.101

param <- c(b0, b1, pt, rt)
xvec <- read.table("xvec.txt")
xvec <- as.matrix(xvec)


ymat  <- gM1nbgen(param, J, R, xvec)
ymatv <- as.numeric(t(ymat))


gM1nb(param, ymat, J, R, xvec)

dM1_nb(ymat, b0, b1, pt, rt, J, R, xvec, TRUE)
dM1_nb_vec(ymatv, b0, b1, pt, rt, J_i, R, xvec, TRUE)

ymatv <- as.numeric(ymat[1, ])

gM1nb(param, ymat[1, , drop = FALSE], J, 1, xvec[1])
dM1_nb(ymat[1, , drop = FALSE], b0, b1, pt, rt, J, 1, xvec[1], TRUE)


set.seed(123)
gM1nbgen(param, J, R, xvec)
set.seed(123)
rM1_nb(1, b0, b1, pt, rt, J, R, xvec)
set.seed(123)
rM1_nb_vec(1, b0, b1, pt, rt, J_i, R, xvec)

set.seed(123)
gM1nbgen(param, J, 1, xvec[1])
set.seed(123)
rM1_nb(1, b0, b1, pt, rt, J, 1, xvec[1])


ymatv_uneven <- rM1_nb_vec(1, b0, b1, pt, rt, J_i_uneven, R, xvec)
dM1_nb_vec(ymatv_uneven, b0, b1, pt, rt, J_i_uneven, R, xvec, TRUE)



mut_i <- b0 + b1 * xvec

nc <- nimbleCode({
   x[1:J] ~ dNmixture_MNB_s(mut = mut, pt = logit(p), rt = rt, J = J)

 })


nmix <- nimbleModel(nc,
                    constants = list(J = J),
                    data = list(x = ymatv[1:J]),
                    inits = list(mut = mut_i[1],
                                 p = p,
                                 rt = rt))
nmix$calculate()

dNmixture_MNB_s(ymatv[1:J], mut_i[1], pt, rt, J, log = 1)
dM0_nb_vec(ymat[1, ], mut_i[1], pt, rt, J, 1, TRUE)



ymat  <- gM1nbgen(param, J, R, xvec)
ymatv <- as.numeric(t(ymat))




nc <- nimbleCode({

   for (i in 1:R) {
 
     mut_i[i] <- b0 + b1 * z[i]
     x[spot1[i]:spot2[i]] ~ dNmixture_MNB_s(mut = mut_i[i], pt = pt, rt = rt, J = J_i[i])
   }
 })


nmix <- nimbleModel(nc,
                    constants = list(J_i = J_i, R = R, spot1 = spot1, spot2 = spot2),
                    data = list(x = ymatv, z = xvec[,1]),
                    inits = list(b0 = b0,
                                 b1 = b1,
                                 pt = pt,
                                 rt = rt))
nmix$calculate()

dM1_nb_vec(ymatv, b0, b1, pt, rt, J_i, R, xvec, TRUE)



nmix <- nimbleModel(nc,
                    constants = list(J_i = J_i_uneven, R = R, spot1 = spot1_uneven, spot2 = spot2_uneven),
                    data = list(x = ymatv_uneven, z = xvec[,1]),
                    inits = list(b0 = b0,
                                 b1 = b1,
                                 pt = pt,
                                 rt = rt))

nmix$calculate()
dM1_nb_vec(ymatv_uneven, b0, b1, pt, rt, J_i_uneven, R, xvec, TRUE)







###################################
#           M2                    #
###################################

source("m2_haines.R")
source("m2_nimble.R")
source("m2_nimble_gen.R")

g0 <- -2
g1 <- 0.25


param <- c(mut, g0, g1, rt)
tvec <- read.table("tvec.txt")
tvec <- as.matrix(tvec)

ymat  <- gM2nbgen(param, J, R, tvec)
ymatv <- as.numeric(t(ymat))

set.seed(123)
gM2nbgen(param, J, R, tvec)
set.seed(123)
rM2_nb(1, mut, g0, g1, rt, J, R, tvec)
set.seed(123)
rM2_nb_vec(1, mut, g0, g1, rt, J_i, R, tvec)

gM2nb(param, ymat, J, R, tvec)
dM2_nb(ymat, mut, g0, g1, rt, J, R, tvec, TRUE)
dM2_nb_vec(ymatv, mut, g0, g1, rt, J_i, R, tvec, TRUE)


ymatv_uneven <- rM2_nb_vec(1, mut, g0, g1, rt, J_i_uneven, R, tvec)
dM2_nb_vec(ymatv_uneven, mut, g0, g1, rt, J_i_uneven, R, tvec, TRUE)



pt_i <- g0 + g1 * tvec[,1]


nc <- nimbleCode({
   x[1:J] ~ dNmixture_MNB_s(mut = mut, pt = pt, rt = rt, J = J)

 })


nmix <- nimbleModel(nc,
                    constants = list(J = J),
                    data = list(x = ymatv[1:J]),
                    inits = list(mut = mut,
                                 pt = pt_i[1],
                                 rt = rt))
nmix$calculate()

dNmixture_MNB_s(ymatv[1:J], mut, pt_i[1], rt, J, log = 1)


ymat  <- gM2nbgen(param, J, R, tvec)
ymatv <- as.numeric(t(ymat))


nc <- nimbleCode({

   for (i in 1:R) {
 
     pt_i[i] <- g0 + g1 * w[i]
     x[spot1[i]:spot2[i]] ~ dNmixture_MNB_s(mut = mut, pt = pt_i[i], rt = rt, J = J_i[i])
   }
 })


nmix <- nimbleModel(nc,
                    constants = list(J_i = J_i, R = R, spot1 = spot1, spot2 = spot2),
                    data = list(x = ymatv, w = tvec[,1]),
                    inits = list(mut = mut,
                                 g0 = g0,
                                 g1 = g1,
                                 rt = rt))
nmix$calculate()

dM2_nb_vec(ymatv, mut, g0, g1, rt, J_i, R, tvec, TRUE)



nmix <- nimbleModel(nc,
                    constants = list(J_i = J_i_uneven, R = R, spot1 = spot1_uneven, spot2 = spot2_uneven),
                    data = list(x = ymatv_uneven, w = tvec[,1]),
                    inits = list(mut = mut,
                                 g0 = g0,
                                 g1 = g1,
                                 rt = rt))

nmix$calculate()
dM2_nb_vec(ymatv_uneven, mut, g0, g1, rt, J_i_uneven, R, tvec, TRUE)



###################################
#           M3                    #
###################################

source("m3_haines.R")
source("m3_nimble.R")
source("m3_nimble_gen.R")
source("m3_nimbleEcology.R")


param <- c(mut, g0, g1, rt)
tvec <- read.table("tvec3.txt")
tvec <- as.matrix(tvec)
tvec_even <- rep(tvec, length.out = max(J_i))[sequence(J_i)]


ymat  <- gM3nbgen(param, J, R, tvec)
ymatv <- as.numeric(t(ymat))


set.seed(123)
gM3nbgen(param, J, R, tvec)
set.seed(123)
rM3_nb(1, mut, g0, g1, rt, J, R, tvec)
set.seed(123)
rM3_nb_vec(1, mut, g0, g1, rt, J_i, R, tvec_even)


gM3nb(param, ymat, J, R, tvec)
dM3_nb(ymat, mut, g0, g1, rt, J, R, tvec, TRUE)
dM3_nb_vec(ymatv, mut, g0, g1, rt, J_i, R, tvec_even, TRUE)

tvec_uneven <- rep(tvec, length.out = max(J_i_uneven))[sequence(J_i_uneven)]
ymatv_uneven <- rM3_nb_vec(1, mut, g0, g1, rt, J_i_uneven, R, tvec_uneven)
dM3_nb_vec(ymatv_uneven, mut, g0, g1, rt, J_i_uneven, R, tvec_uneven, TRUE)




pt_j <- g0 + g1 * tvec[,1]


nc <- nimbleCode({
   x[1:J] ~ dNmixture_MNB_v(mut = mut, pt = pt[1:3], rt = rt, J = J)

 })


nmix <- nimbleModel(nc,
                    constants = list(J = J),
                    data = list(x = ymatv[1:J]),
                    inits = list(mut = mut,
                                 pt = pt_j,
                                 rt = rt))
nmix$calculate()

dNmixture_MNB_v(ymatv[1:J], mut, pt_j, rt, J, log = 1)




nc <- nimbleCode({

   pt_j[1:J_tot] <- g0 + g1 * w[1:J_tot]

   for (i in 1:R) {
 
     x[spot1[i]:spot2[i]] ~ dNmixture_MNB_v(mut = mut, pt = pt_j[spot1[i]:spot2[i]], rt = rt, J = J_i[i])
   }
 })


nmix <- nimbleModel(nc,
                    constants = list(J_i = J_i, R = R, spot1 = spot1, spot2 = spot2, J_tot = J_tot),
                    data = list(x = ymatv, w = tvec_even),
                    inits = list(mut = mut,
                                 g0 = g0,
                                 g1 = g1,
                                 rt = rt))
nmix$calculate()

dM3_nb_vec(ymatv, mut, g0, g1, rt, J_i, R, tvec_even, TRUE)



nmix <- nimbleModel(nc,
                    constants = list(J_i = J_i_uneven, R = R, spot1 = spot1_uneven, spot2 = spot2_uneven, J_tot = J_tot_uneven),
                    data = list(x = ymatv_uneven, w = tvec_uneven),
                    inits = list(mut = mut,
                                 g0 = g0,
                                 g1 = g1,
                                 rt = rt))

nmix$calculate()
dM3_nb_vec(ymatv_uneven, mut, g0, g1, rt, J_i_uneven, R, tvec_uneven, TRUE)




###################################
#           M12                   #
###################################

source("m12_haines.R")
source("m12_nimble.R")


param <- c(b0, b1, g0, g1, rt)
tvec <- read.table("tvec.txt")
tvec <- as.matrix(tvec)


ymat  <- gM12nbgen(param, J, R, xvec, tvec)


gM12nb(param, ymat, J, R, xvec, tvec)
dM12_nb(ymat, b0, b1, g0, g1, rt, J, R, xvec, tvec, TRUE)

set.seed(123)
gM12nbgen(param, J, R, xvec, tvec)
set.seed(123)
rM12_nb(1, b0, b1, g0, g1, rt, J, R, xvec, tvec)

###################################
#           M123                  #
###################################

source("m123_haines.R")
source("m123_nimble.R")
source("m123_nimble_gen.R")

g1 <- 0.125
g2 <- 0.125
param <- c(b0, b1, g0, g1, g2, rt)


tvecR <- read.table("tvec.txt")
tvecR <- as.matrix(tvecR)
tvecR_even <- rep(tvecR, J_i)
tvecJ <- read.table("tvec3.txt")
tvecJ <- as.matrix(tvecJ)
tvecJ_even <- rep(tvecJ, length.out = max(J_i))[sequence(J_i)]


ymat  <- gM123nbgen(param, J, R, xvec, tvecR, tvecJ)
ymatv <- as.numeric(t(ymat))


gM123nb(param, ymat, J, R, xvec, tvecR, tvecJ)
dM123_nb(ymat, b0, b1, g0, g1, g2, rt, J, R, xvec, tvecR, tvecJ, TRUE)
dM123_nb_vec(ymatv, b0, b1, g0, g1, g2, rt, J_i, R, xvec, tvecR_even, tvecJ_even, TRUE)

set.seed(123)
gM123nbgen(param, J, R, xvec, tvecR, tvecJ)
set.seed(123)
rM123_nb(1, b0, b1, g0, g1, g2, rt, J, R, xvec, tvecR, tvecJ)
set.seed(123)
rM123_nb_vec(1, b0, b1, g0, g1, g2, rt, J_i, R, xvec, tvecR_even, tvecJ_even)


tvecJ_uneven <- rep(tvecJ, length.out = max(J_i_uneven))[sequence(J_i_uneven)]
tvecR_uneven <- rep(tvecR, J_i_uneven)
ymatv_uneven <- rM123_nb_vec(1, b0, b1, g0, g1, g2, rt, J_i_uneven, R, xvec, tvecR_uneven, tvecJ_uneven)
dM123_nb_vec(ymatv_uneven, b0, b1, g0, g1, g2, rt, J_i_uneven, R, xvec, tvecR_uneven, tvecJ_uneven, TRUE)







nc <- nimbleCode({

   pt_j[1:J_tot] <- g0 + g1 * wR[1:J_tot] + g2 * wJ[1:J_tot]

   for (i in 1:R) {
     mut_i[i] <- b0 + b1 * z[i] 
     x[spot1[i]:spot2[i]] ~ dNmixture_MNB_v(mut = mut_i[i], pt = pt_j[spot1[i]:spot2[i]], rt = rt, J = J_i[i])
   }
 })


nmix <- nimbleModel(nc,
                    constants = list(J_i = J_i, R = R, spot1 = spot1, spot2 = spot2, J_tot = J_tot),
                    data = list(x = ymatv, z = xvec[,1], wR = tvecR_even, wJ = tvecJ_even),
                    inits = list(b0 = b0,
                                 b1 = b1,
                                 g0 = g0,
                                 g1 = g1,
                                 g2 = g2,
                                 rt = rt))
nmix$calculate()

dM123_nb_vec(ymatv, b0, b1, g0, g1, g2, rt, J_i, R, xvec, tvecR_even, tvecJ_even, TRUE)



nmix <- nimbleModel(nc,
                    constants = list(J_i = J_i_uneven, R = R, spot1 = spot1_uneven, spot2 = spot2_uneven, J_tot = J_tot_uneven),
                    data = list(x = ymatv_uneven, z = xvec[,1], wR = tvecR_uneven, wJ = tvecJ_uneven),
                    inits = list(b0 = b0,
                                 b1 = b1,
                                 g0 = g0,
                                 g1 = g1,
                                 g2 = g2,
                                 rt = rt))

nmix$calculate()
dM123_nb_vec(ymatv_uneven, b0, b1, g0, g1, g2, rt, J_i_uneven, R, xvec, tvecR_uneven, tvecJ_uneven, TRUE)





