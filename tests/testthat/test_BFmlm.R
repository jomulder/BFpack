#several tests if output is the same

# testing coefficients in multivariate normal model
lm1 <- lm(cbind(mpg,cyl,hp) ~ disp + wt, data = mtcars)
BF1 <- BF(lm1)
PHPexplo <- matrix(
  c(0.000,  0.000,  1.000,
  0.373,  0.604,  0.023,
  0.087,  0.908,  0.004,
  0.000,  0.000,  1.000,
  0.000,  0.000,  1.000,
  0.741,  0.178,  0.082,
  0.284,  0.017,  0.700,
  0.007,  0.000,  0.993,
  0.697,  0.239,  0.064),nrow=9,byrow=T)
expect_equivalent(
  round(BF1$PHP_exploratory,3),PHPexplo
)
# tests on same predictor on different DVs
set.seed(123)
BF2 <- BF(x=lm1,hypothesis="disp_on_mpg>disp_on_cyl>disp_on_hp>0; disp_on_mpg=disp_on_cyl=disp_on_hp=0")
test_that("BF.mlm two hypotheses for same IVs correctly evaluated", {
  expect_equivalent(
    round(BF2$PHP_confirmatory,3),c(0.016,0.000,0.984)
)})
# tests on different predictors on same DVs
BF3 <- BF(lm1,hypothesis="disp_on_mpg>wt_on_mpg;disp_on_mpg=wt_on_mpg;disp_on_mpg<wt_on_mpg")
test_that("BF.mlm two hypotheses on same DVs correctly evaluated", {
expect_equivalent(
  round(BF3$PHP_confirmatory,3),c(0.902,0.094,0.005)
)})
# tests on different predictors on different DVs
BF4 <- BF(lm1,hypothesis="disp_on_mpg<wt_on_cyl & disp_on_cyl<wt_on_hp; disp_on_mpg=wt_on_cyl & disp_on_cyl=wt_on_hp")
test_that("BF.mlm two hypotheses different DVs and different IVs correctly evaluated", {
  expect_equivalent(
    round(BF4$PHP_confirmatory,3),c(0.018,0.919,0.062)
)})


# testing correlations in multivariate normal model
lm1 <- lm(cbind(mpg,cyl,hp,qsec) ~ disp + wt, data = mtcars)
set.seed(123)
BF1 <- BF(lm1,parameter="correlation")
PHPexplo <- matrix(
c(0.340,  0.590,  0.070,
0.354,  0.571,  0.074,
0.307,  0.054,  0.639,
0.403,  0.119,  0.478,
0.314,  0.620,  0.066,
0.168,  0.800,  0.032),nrow=6,byrow=T)
test_that("BF.mlm exploratory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    BF1$PHP_exploratory,PHPexplo
)})
# (dummy) hypotheses on the correlations
BF2 <- BF(lm1,parameter="correlation",hypothesis="cyl_with_mpg<hp_with_mpg<hp_with_cyl<0;
   cyl_with_mpg<hp_with_mpg<hp_with_cyl<0")
BFtable <- matrix(
  c(1, 0.0207385,     1, 0.0562496,    1, 2.7123315, 2.7123315, 0.4225317,
    1, 0.0205368,     1, 0.0562507,    1, 2.7390162, 2.7390162, 0.4266887,
    1, 0.9779000,     1, 0.9465000,    1, 0.9678904, 0.9678904, 0.1507796),nrow=3,byrow=T)
test_that("BF.mlm confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    round(BF2$BFtable_confirmatory,7),BFtable
)})


