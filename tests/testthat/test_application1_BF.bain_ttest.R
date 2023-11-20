
# Frequentist one sample t test at the null point mu=5
ttest1 <- t_test(therapeutic,alternative="greater",mu=5)
# one-sided Bayesian one sample t test
BF1 <- BF(ttest1,hypothesis="mu<5",prior.hyp=c(.5,.5))
test_that("bartlett_htest works for a specific one-sided test", {
  expect_equivalent(
  round(ttest1$p.value,3),round(BF1$PHP_confirmatory[1],3)
)})
BF1b <- BF(ttest1,hypothesis="mu=5;mu<5")
test_that("bartlett_htest works the same for confirmatory and explotory setting", {
  expect_equivalent(
    round(BF1b$PHP_exploratory,3),round(BF1b$PHP_confirmatory,3)
  )})

# confirmatory Bayesian one sample t test
BF2 <- BF(ttest1,hypothesis="mu=5;mu>5",prior.hyp=c(.5,.5,0))
#check prior output
test_that("bartlett_htest works for a multiple hypothesis test", {
expect_equivalent(
  BF2$prior,c(.5,.5,0)
)})
#check posterior probabilities
test_that("bartlett_htest works for a multiple hypothesis test", {
expect_equivalent(
  round(BF2$PHP_confirmatory,3),c(0.943,0.057,0.000)
)
})

test_that("t_test prior correct in output", {expect_equivalent(
  BF2$prior.hyp,c(.5,.5,0)
)})
#check posterior probabilities
test_that("t_test PHP correct", {expect_equivalent(
  round(BF2$PHP_confirmatory,3),c(0.943,0.057,0.000)
)})



