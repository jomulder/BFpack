
# test for BF on rma.uni object
library(metafor)

# ### Generate data
tau2 <- 0.05 # True between-study variance

set.seed(123)
vi <- runif(50, min = 0.01, max = 0.2) # Observed sampling variances
yi <- rnorm(50, mean = 0, sd = sqrt(vi+tau2)) # Observed effect sizes

test_that("exploratory metafor test for I^2", {
  skip_on_cran()
  ### Fit a random-effects model to the data
  res <- rma(yi = yi, vi = vi)
  set.seed(123)
  BFmeta1 <- BF(res)
  set.seed(123)
  BFmeta1a <- BF(res,prior.hyp.explo = c(5:7))
  expect_equivalent(
    round(BFmeta1$PHP_exploratory[1,],3),c(0.086,0.013,0.902), tolerance = .05
  )
  expect_equivalent(
    round(BFmeta1$PHP_exploratory[2,],3),c(.800,0.086,0.114), tolerance = .05
  )
  expect_equivalent(
    unname(round(BFmeta1a$PHP_exploratory[2,],3)),
    unname(round(BFmeta1$BFtu_exploratory[2,]*(5:7)/sum(BFmeta1$BFtu_exploratory[2,]*(5:7)),3)),
    tolerance = .05
  )
})

# test fixed effects meta model
test_that("exploratory metafor test for I^2", {
  res <- metafor::rma(yi = yi, vi = vi, method = "EE")
  BFmeta2 <- BF(res)
  expect_equivalent(
    round(BFmeta2$PHP_exploratory,3),c(0.844,0.082,0.074), tolerance = .05
  )
})



