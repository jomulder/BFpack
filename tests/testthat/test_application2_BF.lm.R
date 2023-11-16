aov1 <- aov(price ~ anchor * motivation, data = tvprices)
#names(aov1$coefficients)
BF1 <- BF(aov1)
#check number of tests in exploratory ANOVA
test_that("BF.lm works for an exploratory ANOVA test regarding number of tests", {
expect_equivalent(
  dim(BF1$BFtu_exploratory),c(3,2)
)})
test_that("BF.lm works for a specific exploratory ANOVA test", {
expect_equivalent(
  round(BF1$PHP_exploratory[,1],3),c(0,0,.251)
)})

BF2 <- BF(aov1,hypothesis="anchorrounded = anchorrounded:motivationlow < motivationlow;
          (anchorrounded, anchorrounded:motivationlow) < 0 & motivationlow < 0",
          log = TRUE)
test_that("log(BF.aov) confirmatory ANOVA test", {
  expect_equivalent(
  round(unname(BF2$BFtu_confirmatory),3),c(-40.2, -25.8, 0.043),tol=.1
)})
