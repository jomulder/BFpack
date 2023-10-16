#
# # test for BF on an ergm object
# library(ergm)
# library(Bergm)
# library(BFpack)
#
# seed <- 123
# # florentine data
# data(florentine)
# # example analysis
# ergm_fit <- ergm(flomarriage ~ edges +
#                    kstar(2) +
#                    absdiff("wealth"),
#                  control = control.ergm(seed = seed))
# get_estimates(ergm_fit)
# BFergm.test <- BF(ergm_fit,
#                   hypothesis = "0 = absdiff.wealth > kstar2",
#                   main.iters = 500)
#
# #check if confirmatory test are the same
# test_that("BF.ergm one hypotheses correctly evaluated", {
#   expect_true(
#     all.equal(c(0.2,0.8),
#               unname(BFergm.test$PHP_confirmatory), tolerance = .2)
#   )})
# #check if confirmatory test are the same
# test_that("BF.ergm exploratory is correctly evaluated", {
#   expect_true(
#     all.equal(c(0.1,0,0.9),
#               unname(BFergm.test$PHP_exploratory[2,]), tolerance = .2)
#   )})
#
# # same test with bergm
#
# set.seed(222)
#
# # florentine data
# data(florentine)
# # example analysis
# bergm_fit <- bergm(flomarriage ~ kstar(2) + edges + absdiff("wealth"),
#                  seed = 1,main.iters = 500)
# BFbergm.test <- BF(bergm_fit,hypothesis = "0 = theta3 > theta1",main.iters = 500)
#
# #check if confirmatory test are the same
# test_that("BF.bergm one hypotheses correctly evaluated", {
#   expect_true(
#     all.equal(0.17,
#               unname(BFbergm.test$PHP_confirmatory)[1], tolerance = .2)
#   )})
#
#
