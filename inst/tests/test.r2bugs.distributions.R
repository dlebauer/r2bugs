set.seed(10)
n = 5 # number of tests per distribution
r.distns <- data.frame(distn = rep(c('norm', 'lnorm',
                         'weibull', 'gamma','chisq',
                         'binom', 'nbinom'), each = n),
                       parama = c(runif(n * 2, -10, 10), # norm, lnorm: real
                         runif(n * 3, 0, 3),              # weibull, gamma, chisq: > 0
                         sample(1:(10 * n), size = 2 * n)), #bin, negbin
                       paramb = c(runif(4 * n, 0, 5), #norm, lnorm, weib, gam
                         rep(NA, n), # chisq
                         runif(2 * n, 0, 1)))

bugs.distns               <-  r2bugs.distributions(r.distns)


test_that("r2bugs converts distributions as expected", {
  ## transforms are the same as in function
  expect_true(all(c("weib", "bin", "chisqr", "negbin") %in% bugs.distns$distn))
  norm <- bugs.distns$distn %in% c("norm", "lnorm")
  expect_equal(bugs.distns$paramb[norm],
               1 / r.distns$paramb[norm]^2)
  weib <- bugs.distns$distn == "weib"
  expect_equal(bugs.distns$paramb[weib],
               1 / r.distns$paramb[weib] ^ r.distns$parama[weib])
  bin <- grepl("bin", bugs.distns$distn)
  expect_equal(bugs.distns$paramb[bin],
               r.distns$parama[bin])
  expect_equal(bugs.distns$parama[bin],
               r.distns$paramb[bin])
})

r.distns.recovered.bugs2r <-  bugs2r.distributions(bugs.distns)
r.distns.recovered.r2bugs <-  r2bugs.distributions(bugs.distns, direction = "bugs2r")

test_that("bugs2r.distributions == r2bugs.distributions(..., direction = 'bugs2r')",{
  expect_equal(r.distns.recovered.bugs2r,
               r.distns.recovered.r2bugs)
})

test_that("bugs2r converts distributions back, as expected", {
  r.distns.recovered <- r.distns.recovered.bugs2r
  expect_equal(as.character(r.distns$distn),
               as.character(r.distns.recovered$distn))
  expect_equal(r.distns$parama,
               r.distns.recovered$parama)
  expect_equal(r.distns$paramb[!is.na(r.distns$paramb)],
                 r.distns.recovered$paramb[!is.na(r.distns$paramb)])
})


## Now implement some of these distributions in JAGS

test.r2bugs <- function(r.distn = data.frame(distn = "norm", parama = 0, paramb = 10),
                        n.iter = 100000) {
  bugs.dist <- r2bugs.distributions(r.distn)
  Y.BUGS <- bugs.rdist(bugs.dist, n.iter = n.iter)
  if(!grepl("chisq", r.distn$distn)){
    Y.R    <- do.call(paste("r", r.distn$distn, sep = ""), list(n.iter/4, r.distn$parama, r.distn$paramb))
  } else {
    Y.R    <- do.call(paste("r", r.distn$distn, sep = ""), list(n.iter/4, r.distn$parama))
  }
  median.test <- (median(Y.BUGS) - median(Y.R)) / (median(c(Y.BUGS, Y.R)))
  var.test <- (var(Y.BUGS) - var(Y.R)) / (var(c(Y.BUGS, Y.R)))
  return(list(median.test = median.test, var.test = var.test, Y.BUGS = Y.BUGS, Y.R = Y.R))
}

test_that("JAGS creates reasonable parameters from these distributions",{
  set.seed(0)
  for(i in 31:nrow(r.distns)) {
    d <- r.distns[i,]
    print(d)
    print(r2bugs.distributions(d))
    ans <- test.r2bugs(r.distn = d)
    print(mean(ans$Y.BUGS))
    print(mean(ans$Y.R))
    expect_true(abs(ans$median.test) < 0.5)
    expect_true(abs(ans$var.test)  < 0.5 | var(ans$Y.R) > 1e10)
  }
})
