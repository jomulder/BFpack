
# test for BF on rma.uni object
library(metafor)

# ### Generate data
tau2 <- 0.05 # True between-study variance

set.seed(123)
vi <- runif(50, min = 0.01, max = 0.2) # Observed sampling variances
yi <- rnorm(50, mean = 0, sd = sqrt(vi+tau2)) # Observed effect sizes

test_that("exploratory metafor random effects model", {
  skip_on_cran()
  ### Fit a random-effects model to the data
  res <- rma(yi = yi, vi = vi)
  set.seed(123)
  BFmeta1 <- BF(res,BF.type="stand.effect",cov.prob=.9)
  set.seed(123)
  BFmeta1a <- BF(res,BF.type="stand.effect",prior.hyp.explo = c(5:7))
  set.seed(123)
  BFmeta1b <- BF(res,BF.type=prior("norm", c(mean = 0, sd = .5)))
  expect_equivalent(
    round(BFmeta1$PHP_exploratory[1,],3),c(.950,.020,.039), tolerance = .05
  )
  expect_equivalent(
    round(BFmeta1b$PHP_exploratory[1,],3),c(.906,.055,.039), tolerance = .05
  )
  expect_equivalent(
    unname(round(BFmeta1a$PHP_exploratory[2,],3)),
    unname(round(BFmeta1$BFtu_exploratory[2,]*(5:7)/sum(BFmeta1$BFtu_exploratory[2,]*(5:7)),3)),
    tolerance = .05
  )
})

# test fixed effects meta model
test_that("exploratory metafor fixed effects model", {
  set.seed(123)
  res <- metafor::rma(yi = yi, vi = vi, method = "EE")
  BFmeta2 <- BF(res,BF.type="stand.effect")
  expect_equivalent(
    round(BFmeta2$PHP_exploratory,3),c(.964,0.019,0.017), tolerance = .05
  )
  set.seed(123)
  res$ni <- rep(10,50)
  BFmeta2b <- BF(res,BF.type="unit.info")
  expect_equivalent(
    round(BFmeta2b$PHP_exploratory,3),c(.957,0.023,0.022), tolerance = .05
  )
})



