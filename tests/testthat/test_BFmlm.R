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
BF2 <- BF(x=lm1,hypothesis="disp_on_mpg>disp_on_cyl>disp_on_hp>0;disp_on_mpg=disp_on_cyl=disp_on_hp=0")
expect_equivalent(
  round(BF2$PHP_confirmatory,3),c(0.016,0.000,0.984)
)
# tests on different predictors on same DVs
BF3 <- BF(lm1,hypothesis="disp_on_mpg>wt_on_mpg;disp_on_mpg=wt_on_mpg;disp_on_mpg<wt_on_mpg")
expect_equivalent(
  round(BF3$PHP_confirmatory,3),c(0.902,0.094,0.005)
)
# tests on different predictors on different DVs
BF4 <- BF(lm1,hypothesis="disp_on_mpg<wt_on_cyl & disp_on_cyl<wt_on_hp; disp_on_mpg=wt_on_cyl & disp_on_cyl=wt_on_hp")
expect_equivalent(
  round(BF4$PHP_confirmatory,3),c(0.018,0.919,0.062)
)


# testing correlations in multivariate normal model
lm1 <- lm(cbind(mpg,cyl,hp,qsec) ~ disp + wt, data = mtcars)
set.seed(123)
BF1 <- BF(lm1,parameter="correlation")
PHPexplo <- matrix(
c(0.334,  0.600,  0.067,
0.369,  0.551,  0.080,
0.309,  0.054,  0.638,
0.411,  0.120,  0.469,
0.336,  0.587,  0.078,
0.161,  0.809,  0.030),nrow=6,byrow=T)
expect_equivalent(
  BF1$PHP_exploratory,PHPexplo
)
# (dummy) hypotheses on the correlations
BF2 <- BF(lm1,parameter="correlation",hypothesis="cyl_with_mpg<hp_with_mpg<hp_with_cyl<0;
   cyl_with_mpg<hp_with_mpg<hp_with_cyl<0")
BFtable <- matrix(
c(1, 0.0208456,     1, 0.0522040,    1, 2.5043164, 2.5043164, 0.4189960,
1, 0.0208492,     1, 0.0521888,    1, 2.5031516, 2.5031516, 0.4188012,
1, 0.9796000,     1, 0.9497000,    1, 0.9694773, 0.9694773, 0.1622028),nrow=3,byrow=T)
expect_equivalent(
  round(BF2$BFtable_confirmatory,7),BFtable
)


