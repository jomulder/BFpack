
# test for BF on rma.uni object
library(metafor)

# ### Generate data
tau2 <- 0.05 # True between-study variance

set.seed(123)
vi <- runif(50, min = 0.01, max = 0.2) # Observed sampling variances
yi <- rnorm(50, mean = 0, sd = sqrt(vi+tau2)) # Observed effect sizes

### Fit a random-effects model to the data
res <- rma(yi = yi, vi = vi)
BFmeta1 <- BF(res)

test_that("exploratory metafor test for I^2", {
  expect_equivalent(
    round(BFmeta1$PHP_exploratory[1,],3),c(0.086,0.013,0.902), tolerance = .05
  )})

test_that("exploratory metafor test for delta", {
  expect_equivalent(
    round(BFmeta1$PHP_exploratory[2,],3),c(.800,0.086,0.114), tolerance = .05
  )})


# test fixed effects meta model
res <- metafor::rma(yi = yi, vi = vi, method = "EE")
BFmeta2 <- BF(res)
test_that("exploratory metafor test for I^2", {
  expect_equivalent(
    round(BFmeta2$PHP_exploratory,3),c(0.844,0.082,0.074), tolerance = .05
  )})

# 1 + 1 = 3



