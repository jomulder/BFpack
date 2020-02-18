#several tests if output is the same
set.seed(36)
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
  0.697,  0.239,  0.064,
  0.346,  0.592,  0.062,
  0.297,  0.659,  0.044,
  0.504,  0.230,  0.266),nrow=12,byrow=T)
test_that("BF.mlm exploratory correctly evaluated for given data sets", {
  expect_equivalent(
    round(BF1$PHP_exploratory,3),PHPexplo, tolerance = .03
  )})

# tests on same predictor on different DVs
set.seed(123)
BF2 <- BF(x=lm1,hypothesis="disp_on_mpg>disp_on_cyl>disp_on_hp>0; disp_on_mpg=disp_on_cyl=disp_on_hp=0")
test_that("BF.mlm two hypotheses for same IVs correctly evaluated", {
  expect_equivalent(
    round(BF2$PHP_confirmatory,3),c(0.016,0.000,0.984), tolerance = .015
)})
# tests on different predictors on same DVs
set.seed(574)
BF3 <- BF(lm1,hypothesis="disp_on_mpg>wt_on_mpg;disp_on_mpg=wt_on_mpg;disp_on_mpg<wt_on_mpg")
test_that("BF.mlm two hypotheses on same DVs correctly evaluated", {
expect_equivalent(
  round(BF3$PHP_confirmatory,3),c(0.902,0.094,0.005)
)})
# tests on different predictors on different DVs
set.seed(4768)
BF4 <- BF(lm1,hypothesis="disp_on_mpg<wt_on_cyl & disp_on_cyl<wt_on_hp; disp_on_mpg=wt_on_cyl & disp_on_cyl=wt_on_hp")
test_that("BF.mlm two hypotheses different DVs and different IVs correctly evaluated", {
  expect_equivalent(
    round(BF4$PHP_confirmatory,3),c(0.018,0.919,0.062), tolerance = .005
)})


# testing correlations in multivariate normal model
lm1 <- lm(cbind(mpg,cyl,hp,qsec) ~ disp + wt, data = mtcars)
set.seed(123)
BF1 <- BF(lm1)
PHPexplo <- matrix(
c(0.340,  0.590,  0.070,
0.354,  0.571,  0.074,
0.307,  0.054,  0.639,
0.403,  0.119,  0.478,
0.314,  0.620,  0.066,
0.168,  0.800,  0.032),nrow=6,byrow=T)
test_that("BF.mlm exploratory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    BF1$PHP_exploratory[-(1:12),],PHPexplo, tolerance = .2
)})
# (dummy) hypotheses on the correlations
set.seed(564574)
BF2 <- BF(lm1,hypothesis="cyl_with_mpg<hp_with_mpg<hp_with_cyl<0;
   cyl_with_mpg<hp_with_mpg<hp_with_cyl<0")
BFtable <- structure(c(1, 1, 1, 0.0204763291181508, 0.0208215890818035,
                       0.9815, 1, 1, 1, 0.0546186204667969, 0.0545440643456688, 0.9447,
                       1, 1, 1, 2.66740293885887, 2.6195918155611, 0.962506367804381,
                       2.3851, 2.4045, 0.962506367804381, 0.426818539062759,
                       0.419168148677558, 0.154013312259683), .Dim = c(3L, 8L), .Dimnames = list(
                         c("H1", "H2", "H3"), c("comp_E", "comp_O", "fit_E", "fit_O",
                                                "BF_E", "BF_O", "BF", "PHP")))

test_that("BF.mlm confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    as.vector(unname(BF2$BFtable_confirmatory))[-c(16, 17, 19, 20, 22, 24)], as.vector(unname(BFtable))[-c(16, 17, 19, 20, 22, 24)], tolerance = .01
    )
  expect_equivalent(
    as.vector(unname(BF2$BFtable_confirmatory))[c(22, 24)], as.vector(unname(BFtable))[c(22, 24)], tolerance = .02
  )
  expect_equivalent(
    as.vector(unname(BF2$BFtable_confirmatory))[c(16, 17, 19, 20)], as.vector(unname(BFtable))[c(16, 17, 19, 20)], tolerance = .4
  )
  })


