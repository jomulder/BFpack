#several tests if output is the same
set.seed(36)
# testing coefficients in multivariate normal model
lm1 <- lm(cbind(mpg,cyl,hp) ~ disp + wt, data = mtcars)
BF1 <- BF(lm1, BF.type = "AFBF",prior.hyp.explo=c(1,1,1))
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
test_that("BF.mlm exploratory correctly evaluated for given data sets", {
  expect_equivalent(
    round(BF1$PHP_exploratory,3),PHPexplo, tolerance = .03
  )})

# tests on same predictor on different DVs
set.seed(123)
BF2 <- BF(x=lm1,hypothesis="disp_on_mpg>disp_on_cyl>disp_on_hp>0; disp_on_mpg=disp_on_cyl=disp_on_hp=0")
test_that("BF.mlm two hypotheses for same IVs correctly evaluated", {
  expect_equivalent(
    round(BF2$PHP_confirmatory,3),c(0.010,0.000,0.988), tolerance = .015
)})
# tests on different predictors on same DVs
set.seed(574)
BF3 <- BF(lm1,hypothesis="disp_on_mpg>wt_on_mpg;disp_on_mpg=wt_on_mpg;disp_on_mpg<wt_on_mpg",
          BF.type="AFBF")
test_that("BF.mlm two hypotheses on same DVs correctly evaluated", {
expect_equivalent(
  round(BF3$PHP_confirmatory,3),c(0.902,0.094,0.005)
)})
BF3.log <- BF(lm1,hypothesis="disp_on_mpg>wt_on_mpg;disp_on_mpg=wt_on_mpg;disp_on_mpg<wt_on_mpg",log=TRUE,
              BF.type="AFBF")
test_that("BF.mlm two hypotheses on same DVs correctly evaluated", {
  expect_equivalent(
    round(BF3.log$BFtu_exploratory,3),round(log(BF3$BFtu_exploratory),3),tol=.01
  )})
test_that("BF.mlm two hypotheses different DVs and different IVs correctly evaluated", {
  skip_on_cran()
  # tests on different predictors on different DVs
  set.seed(4768)
  BF4 <- BF(lm1,hypothesis="disp_on_mpg<wt_on_cyl & disp_on_cyl<wt_on_hp; disp_on_mpg=wt_on_cyl & disp_on_cyl=wt_on_hp",
            log = TRUE, BF.type="AFBF")
  expect_equivalent(
    round(BF4$PHP_confirmatory,3),c(0.018,0.919,0.062), tolerance = .005
  )
  expect_equivalent(
    BF4$PHP_confirmatory[1]/BF4$PHP_confirmatory[2],exp(BF4$BFtu_confirmatory[1] -
            BF4$BFtu_confirmatory[2]), tolerance = .01
  )
})


