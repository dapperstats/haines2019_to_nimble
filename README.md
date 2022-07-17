# Converting code from Haines 2019 to nimble 

[Haines 2019](https://onlinelibrary.wiley.com/doi/10.1111/biom.13147) provides likelihood (density) calculation and data generating (randomization) functions in R for multinomial N-mixture models.
Here, I translate those models into [nimble](https://r-nimble.org) code for use alongside the binomial N-mixture models available in [nimbleEcology](https://github.com/nimble-dev/nimbleEcology). 

Haines' work extends the fast-algorithm approach that nimbleEcology uses with binomial N-mixture models by use of recursion to avoid infinite sums.
