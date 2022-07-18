library(nimble)

source("m0_haines.R")
source("m0_nimble.R")
source("m0_nimble_gen.R")

source("m1_haines.R")

# work space for playing out model dev

R  <- 20
J  <- 3
J_i <- rep(J, R)
set.seed(125)

pt  <- 0.615
mut <- 28.048
rt  <- 11.415

param <- c(mut, pt, rt)
ymat  <- gM0nbgen(param)
ymatv <- as.numeric(t(ymat))


  yrow      <- as.matrix(rowSums(ymat))
  ycol      <- as.matrix(colSums(ymat))
  ytot      <- sum(sum(ymat))
  jvec      <- seq(0, J - 1)
  ysumj     <- sum(ycol * jvec)
  ylogfact  <- sum(sum(log(factorial(ymat))))
  yrlogfact <- sum(log(factorial(yrow)))

param0 <- c(log(mut), logit(pt), log(rt))
gM0nb(param0)


dM0_nb(ymat, mut, pt, rt, J, R, TRUE)
dM0_nb_vec(ymatv, mut, pt, rt, J_i, R, TRUE)



pt  <- 0.615
b0t <- 3.077
b1t <- 0.101
rt  <- 4
param <- c(b0t, b1t, pt, rt)
xvec <- read.table("xvec50.txt")
xvec <- as.matrix(xvec)

R <- 50
J <- 3
ymat  <- gM1nbgen(param)

  yrow      <- rowSums(ymat)
  ycol      <- colSums(ymat)
  ytot      <- sum(sum(ymat))
  jvec      <- seq(0, J - 1)
  ysumj     <- sum(ycol * jvec)
  ylogfact  <- sum(sum(log(factorial(ymat))))
  yrlogfact <- sum(log(factorial(yrow)))

param0 <- c(b0t, b1t, logit(pt), log(rt))
gM1nb(param0)


param0p <- c(logit(pt))
gM1pc(param0p)

