
# test for BF on an ergm object

#check if confirmatory test are the same for ergm object
test_that("BF.ergm tests", {
  skip_on_cran()
  library(ergm)
  seed <- 123
  # florentine data
  data(florentine)
  ergm_fit <- ergm(flomarriage ~ edges +
                     kstar(2) +
                     absdiff("wealth"),
                   control = control.ergm(seed = seed))
  get_estimates(ergm_fit)
  seed <- 123
  BFergm.test <- BF(ergm_fit,
                    hypothesis = "0 = absdiff.wealth > kstar2",
                    main.iters = 500, prior.hyp.explo = c(1,1,1))
  expect_true(
    all.equal(c(0.185,0.815),
              unname(BFergm.test$PHP_confirmatory), tolerance = .2)
  )
  expect_true(
    all.equal(c(0.1,0,0.9),
              unname(BFergm.test$PHP_exploratory[2,]), tolerance = .2)
  )
  seed <- 123
  BFergm.test2 <- BF(ergm_fit,
                     hypothesis = "0 = absdiff.wealth > kstar2",
                     main.iters = 500, prior.hyp.explo = 1:3)
  expect_equivalent(
    unname(BFergm.test2$PHP_exploratory[1,]),
    unname(BFergm.test2$BFtu_exploratory[1,]*(1:3)/sum(BFergm.test2$BFtu_exploratory[1,]*(1:3))),
    tol=.05
  )
})

# same test with bergm
#check if confirmatory test are the same
test_that("BF.bergm one hypotheses correctly evaluated", {
  skip_on_cran()
  # example analysis
  library(Bergm)
  library(ergm)
  data(florentine)
  set.seed(222)
  bergm_fit <- bergm(flomarriage ~ kstar(2) + edges + absdiff("wealth"),
                     seed = 1,main.iters = 500)
  BFbergm.test <- BF(bergm_fit,hypothesis = "0 = theta3 > theta1",main.iters = 500)
  expect_true(
    all.equal(0.17,
              unname(BFbergm.test$PHP_confirmatory)[1], tolerance = .2)
  )
})


