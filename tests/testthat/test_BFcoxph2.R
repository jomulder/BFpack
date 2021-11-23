library(survival)
# Create the simplest test data set
test1 <- list(time=c(4,3,1,1,2,2,3),
              status=c(1,1,1,0,1,1,0),
              x=c(0,2,1,1,1,0,0),
              sex=c(0,0,0,0,1,1,1))
# Fit a stratified model
fit <- coxph(Surv(time, status) ~ x + strata(sex), test1)

out <- BF(fit, "x > .5")

test_that("BF.coxph returns correct results,", {
  expect_true(out$PHP_confirmatory[1] > .5 & out$PHP_confirmatory[2] < .5)
})


# Create a simple data set for a time-dependent model
test2 <- list(start=c(1,2,5,2,1,7,3,4,8,8),
              stop=c(2,3,6,7,8,9,9,9,14,17),
              event=c(1,1,1,1,1,1,1,0,0,0),
              pred=c(1,0,0,1,0,1,1,1,0,0))
fit <- coxph(Surv(start, stop, event) ~ pred, test2)

out <- BF(fit, hypothesis = "pred > .5; pred = .5")
test_that("BF.coxph returns correct results,", {
  expect_true(out$PHP_confirmatory[1] < .5 & out$PHP_confirmatory[2] > .5)
})

out <- BF(fit, hypothesis = "pred > .5; pred = .5",complement = FALSE)
test_that("BF.coxph returns correct results,", {
  expect_equivalent(
    out$PHP_confirmatory, c(0.194,0.806), tol = .1
  )
})
