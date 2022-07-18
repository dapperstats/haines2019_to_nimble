library(nimble)

# M0

source("m0_haines.R")
source("m0_nimble.R")
source("m0_nimble_gen.R")
source("m0_nimbleEcology.R")

R  <- 1
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

ymat  <- matrix(gM0nbgen(param),nrow = 1, ncol = 3)
ymatv <- as.numeric(t(ymat))


  yrow      <- as.matrix(rowSums(ymat))
  ycol      <- as.matrix(colSums(ymat))
  ytot      <- sum(sum(ymat))
  jvec      <- seq(0, J - 1)
  ysumj     <- sum(ycol * jvec)
  ylogfact  <- sum(sum(log(factorial(ymat))))
  yrlogfact <- sum(log(factorial(yrow)))

gM0nb(param)

dM0_nb(ymat, mu, p, r, J, R, TRUE)
dM0_nb_vec(ymatv, mu, p, r, J_i, R, TRUE)


dNmixture_MNB_s(ymat[1, ], mu, p, r, J, TRUE)
dM0_nb(ymat[1, , drop = FALSE], mu, p, r, J, 1, TRUE)


# M1

source("m1_haines.R")
source("m1_nimble.R")


p   <- 0.61
b0  <- 21.693
b1  <- 0.101
r   <- 4
pt  <- logit(p)
rt  <- log(r)
b0t <- log(b0)
b1t <- b1


param <- c(b0t, b1t, pt, rt)
xvec <- read.table("xvec50.txt")
xvec <- as.matrix(xvec)[1]

R <- 1
J <- 3
ymat  <- matrix(gM1nbgen(param), ncol = 3)

  yrow      <- rowSums(ymat)
  ycol      <- colSums(ymat)
  ytot      <- sum(sum(ymat))
  jvec      <- seq(0, J - 1)
  ysumj     <- sum(ycol * jvec)
  ylogfact  <- sum(sum(log(factorial(ymat))))
  yrlogfact <- sum(log(factorial(yrow)))

param0 <- c(b0t, b1t, pt, rt)
gM1nb(param0)

dM1_nb(ymat, b0, b1, p, r, J, R, xvec,TRUE)

ymatv <- as.numeric(ymat)

dNmixture_MNB_sitecovar_s(ymatv, b0, b1, p, r, J, xvec,TRUE)





