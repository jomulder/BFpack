# Test class coxph with Gaussian_estimator
test2 <- list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
              stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
              event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
              x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) )
coxph1 <- survival::coxph(survival::Surv(start, stop, event) ~ x, test2)
BF1 <- BF(coxph1,prior.hyp.explo = c(1,1,1))
test_that("BF.coxph exploratory hypotheses correctly evaluated", {
  expect_true(
    all.equal(round(c(unname(BF1$PHP_exploratory)),3),c(0.569,0.220,0.211))
)})
BF1log <- BF(coxph1,log = TRUE, prior.hyp.explo = c(1,1,1))
test_that("BF.coxph exploratory hypotheses correctly evaluated", {
  expect_true(
    all.equal(exp(c(BF1log$BFtu_exploratory)),c(BF1$BFtu_exploratory))
  )})
test_that("BF.coxph exploratory hypotheses correctly evaluated", {
  expect_true(
    all.equal(c(unname(BF1$PHP_exploratory)),c(unname(BF1log$PHP_exploratory)))
  )})
