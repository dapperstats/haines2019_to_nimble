# Set up constants and initial values for defining the model
J <- 3
dat <- ymat[1,]


# Define code for a nimbleModel
 nc <- nimbleCode({
   x[1:J] ~ dNmixture_MNB_s(mut = log(mu), pt = logit(p),
                            rt = log(r), J = J)

   mu ~ dunif(0, 1000)
   r  ~ dunif(0, 10)
   p  ~ dunif(0, 1)

 })

nmix <- nimbleModel(nc,
                    data = list(x = dat),
                    constants = list(J = J),
                    inits = list(mu = mu,
                                 p = p,
                                 r = r))

nmix$calculate()


dNmixture_s(dat, lambda = mu, prob = p,
                        Nmin = -1, Nmax = -1, len = J)



dNmixture_MNB_s(ymat[1, ], mut, pt, rt, J)