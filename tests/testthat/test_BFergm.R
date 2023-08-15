
# test for BF on an ergm object
library(ergm)
library(Bergm)

seed <- 123
# florentine data
data(florentine)
# example analysis
ergm_fit <- ergm(flomarriage ~ edges + kstar(2) + absdiff("wealth"),
                 control = control.ergm(seed = seed))
#check get.estimates
get_estimates(ergm_fit)

BFergm.test <- BF(ergm_fit,hypothesis = "0 = absdiff.wealth > kstar2")

#check if confirmatory test are the same
test_that("BF.ergm one hypotheses correctly evaluated", {
  expect_true(
    all.equal(c(0.235,0.765),
              unname(BFergm.test$PHP_confirmatory), tolerance = .005)
  )})
#check if confirmatory test are the same
test_that("BF.ergm exploratory is correctly evaluated", {
  expect_true(
    all.equal(c(0.123,0.011,0.866),
              unname(BFergm.test$PHP_exploratory[2,]), tolerance = .005)
  )})

# same test with bergm

set.seed(222)

# florentine data
data(florentine)
# example analysis
bergm_fit <- bergm(flomarriage ~ edges + kstar(2) + absdiff("wealth"),
                 seed = 1)
summary(bergm_fit)

#check get.estimates
get_estimates(bergm_fit)

BFbergm.test <- BF(bergm_fit,hypothesis = "0 = theta3 > theta2")

#check if confirmatory test are the same
test_that("BF.bergm one hypotheses correctly evaluated", {
  expect_true(
    all.equal(c(0.1,0.9),
              unname(BFbergm.test$PHP_confirmatory), tolerance = .2)
  )})


