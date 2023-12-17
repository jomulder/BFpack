
test_that("for mlm object multiple hypotheses with equality and order constraints can be computed", {
  skip_on_cran()
  fmri.lm <- lm(cbind(Superficial,Middle,Deep) ~ Face + Vehicle, data=fmri)
  constraints.fmri <- "Face_on_Deep = Face_on_Superficial = Face_on_Middle < 0 <
    Vehicle_on_Deep = Vehicle_on_Superficial = Vehicle_on_Middle;
    Face_on_Deep < Face_on_Superficial = Face_on_Middle < 0 < Vehicle_on_Deep =
    Vehicle_on_Superficial = Vehicle_on_Middle"
  set.seed(123)
  BF_fmri1 <- BF(fmri.lm, hypothesis = constraints.fmri)
  #check results
  matrixmeasures <- matrix(c(0.2741948, 0.19630283, 0.7186756, 0.9998622,  2.621041,  5.093468,  13.35019, 0.023006434,
                             0.2967301, 0.02525591, 5.1451176, 0.8243133, 17.339388, 32.638428, 565.93036, 0.975270261,
                             1.0000000, 1.00000000, 1.0000000, 1.0000000,  1.000000,  1.000000, 1.00000, 0.001723304),byrow=T,nrow=3)
  expect_equivalent(
    BF_fmri1$BFtable_confirmatory, matrixmeasures, tolerance = .01
  )
})

#no seed required
fmri.lm <- lm(cbind(Superficial,Middle,Deep) ~ Face + Vehicle, data=fmri)
constraints.fmri <- "Face_on_Deep = Face_on_Superficial = Face_on_Middle < 0;
    Face_on_Deep < Face_on_Superficial = Face_on_Middle = 0"
BF_fmri2 <- BF(x = fmri.lm, hypothesis = constraints.fmri)
test_that("for mlm object multiple hypotheses for same IV", {
  expect_equivalent(
    BF_fmri2$PHP_confirmatory,c(0.18941726,0.72354824,0.08703451)
  )})
BF_fmri2b <- BF(fmri.lm, hypothesis = constraints.fmri, log=TRUE)
test_that("for mlm object multiple hypotheses for same IV with log", {
  expect_equivalent(
    exp(BF_fmri2b$BFtu_confirmatory),BF_fmri2$BFtu_confirmatory
)})

#no seed required
constraints.fmri <- "Face_on_Deep = Vehicle_on_Deep < 0;
    Face_on_Middle = Vehicle_on_Middle < 0"
BF_fmri3 <- BF(x = fmri.lm, hypothesis = constraints.fmri)
test_that("for mlm object multiple hypotheses on same DV", {
  expect_equivalent(
    BF_fmri3$PHP_confirmatory,c(0.05595577,0.38093401,0.56311022 )
  )})

# test multivariate t test
# the hypothesis argument does not allow that parameters start with a bracket (because it is also used
# for other reasons), so the standard intercept, i.e., (intercept) is problematic to use.
intercept <- rep(1,nrow(fmri))
lm1 <- lm(cbind(Superficial,Middle,Deep) ~ -1 + intercept, data=fmri)
# test the three means jointly against the null vector (1,1,1)
BF1 <- BF(lm1,hypothesis="intercept_on_Superficial=1 & intercept_on_Middle=1 & intercept_on_Deep=1")
mvt_test(fmri)



