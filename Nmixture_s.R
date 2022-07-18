dNmixture_s <- nimbleFunction(
    run = function(x = double(1),
                   lambda = double(),
                   prob = double(),
                   Nmin = double(0, default = -1),
                   Nmax = double(0, default = -1),
                   len = double(),
                   log = integer(0, default = 0)) {
  if (length(x) != len) stop("in dNmixture_s, len must equal length(x).")

  # Lambda cannot be negative
  if (lambda < 0) {
    if (log) return(-Inf)
    else return(0)
  }

  ## For each x, the conditional distribution of (N - x | x) is pois(lambda * (1-p))
  ## We determine the lowest N and highest N at extreme quantiles and sum over those.
  if (Nmin == -1) {
    Nmin <- min(x + qpois(0.00001, lambda * (1 - prob)))
  }
  if (Nmax == -1) {
    Nmax <- max(x + qpois(0.99999, lambda * (1 - prob)))
  }
  Nmin <- max( max(x), Nmin ) ## set Nmin to at least the largest x

  logProb <- -Inf

  if (Nmax > Nmin) {
    numN <- Nmax - Nmin + 1 - 1  ## remember: +1 for the count, but -1 because the summation should run from N = maxN to N = minN + 1
    prods <- rep(0, numN)
    for (i in (Nmin + 1):Nmax) {
      prods[i - Nmin] <- prod(i/(i - x)) / i
    }

    ff <- log(lambda) + log(1-prob)*len + log(prods)
    log_fac <- nimNmixPois_logFac(numN, ff)
    logProb <- dpois(Nmin, lambda, log = TRUE) + sum(dbinom(x, size = Nmin, prob = prob, log = TRUE)) + log_fac
  }
  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())
})