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
matrix(gM0nbgen(param, J, R), nrow = R, ncol = J)


###################################
#           M1                    #
###################################

source("m1_haines.R")
source("m1_nimble.R")
source("m1_nimbleEcology.R")


b0 <- 3.077
b1 <- 0.101

param <- c(b0, b1, pt, rt)
xvec <- read.table("xvec50.txt")
xvec <- as.matrix(xvec)


ymat  <- gM1nbgen(param, J, R, xvec)


gM1nb(param, ymat, J, R, xvec)

dM1_nb(ymat, b0, b1, pt, rt, J, R, xvec, TRUE)

ymatv <- as.numeric(ymat[1, ])

gM1nb(param, ymat[1, , drop = FALSE], J, 1, xvec[1])
dM1_nb(ymat[1, , drop = FALSE], b0, b1, pt, rt, J, 1, xvec[1], TRUE)
dNmixture_MNB_sitecovar_s(ymatv, b0, b1, pt, rt, J, xvec[1], TRUE)


set.seed(123)
gM1nbgen(param, J, R, xvec)
set.seed(123)
rM1_nb(1, b0, b1, pt, rt, J, R, xvec)

set.seed(123)
gM1nbgen(param, J, 1, xvec[1])
set.seed(123)
rM1_nb(1, b0, b1, pt, rt, J, 1, xvec[1])
set.seed(123)
rNmixture_MNB_sitecovar_s(1, b0, b1, pt, rt, J, xvec[1])




###################################
#           M2                    #
###################################

source("m2_haines.R")

g0 <- -2
g1 <- 0.25


param <- c(mut, g0, g1, rt)
tvec <- read.table("tvec.txt")
tvec <- as.matrix(tvec)

ymat  <- gM2nbgen(param, J, R, tvec)


set.seed(123)
gM2nbgen(param, J, R, tvec)
set.seed(123)
rM2_nb(1, mut, g0, g1, rt, J, R, tvec)

gM2nb(param, ymat, J, R, tvec)
dM2_nb(ymat, mut, g0, g1, rt, J, R, tvec, TRUE)




