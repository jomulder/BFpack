
# test for BF on an ergm object
library(ergm)

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
    all.equal(c(0.194,0.806),
              unname(BFergm.test$PHP_confirmatory), tolerance = .005)
  )})
#check if confirmatory test are the same
test_that("BF.ergm exploratory is correctly evaluated", {
  expect_true(
    all.equal(c(0.413,0.498,0.088),
              unname(BFergm.test$PHP_exploratory[2,]), tolerance = .005)
  )})

