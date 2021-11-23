survreg1 <- survival::survreg(survival::Surv(futime, fustat) ~ ecog.ps + rx, survival::ovarian,
                              dist='weibull', scale=1)
BF1 <- BF(survreg1)
test_that("BF.survreg exploratory hypotheses correctly evaluated", {
  expect_equivalent(
    round(BF1$BFtu_exploratory[,1],5),c(0.00000,3.88378,3.12164)
  )})
