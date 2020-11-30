fit <- glm(sent ~ ztrust + zfWHR + zAfro + glasses + attract + maturity +
             tattoos, family = binomial(), data = wilson)
set.seed(123)
BF_glm <- BF(fit, hypothesis="ztrust > zfWHR > zAfro > 0;
             ztrust > 0 & zfWHR = zAfro = 0")
#check results exploratory test
test_that("glm two hypotheses correctly evaluated via exploratory test", {
expect_equivalent(
  BF_glm$PHP_exploratory[,3],c(0.133,1.000,0.999,0.004,0.278,0.029,0.011,0.202)
)})
#check results conformatory test
test_that("glm two hypotheses correctly evaluated via confirmatory test", {
  expect_equivalent(
  BF_glm$PHP_confirmatory,c(0.142, 0.002, 0.856), tolerance = .02
)})
#check equal prior probabilities
test_that("glm use correct prior probabilities", {
expect_true(
  all.equal(unname(BF_glm$prior.hyp),rep(1,3)/3)
)})

#check that complement is not added the complete parameter space is covered by hypotheses
BF2 <- BF(fit, hypothesis="ztrust = zfWHR;
             ztrust > zfWHR; ztrust < zfWHR")
test_that("glm does not add a complement hypothesis as the test is exhaustive", {
expect_true(
  length(BF2$PHP_confirmatory)==3
)})

