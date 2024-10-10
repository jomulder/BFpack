
test_that("for mlm object multiple hypotheses with equality and order constraints can be computed", {
  skip_on_cran()
  fmri.lm <- lm(cbind(Superficial,Middle,Deep) ~ Face + Vehicle, data=fmri)
  constraints.fmri <- "Face_on_Deep = Face_on_Superficial = Face_on_Middle < 0 <
    Vehicle_on_Deep = Vehicle_on_Superficial = Vehicle_on_Middle;
    Face_on_Deep < Face_on_Superficial = Face_on_Middle < 0 < Vehicle_on_Deep =
    Vehicle_on_Superficial = Vehicle_on_Middle"
  set.seed(123)
  BF_fmri1 <- BF(fmri.lm, hypothesis = constraints.fmri, BF.type = "AFBF")
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
BF_fmri2 <- BF(x = fmri.lm, hypothesis = constraints.fmri, BF.type = "AFBF")
test_that("for mlm object multiple hypotheses for same IV", {
  expect_equivalent(
    BF_fmri2$PHP_confirmatory,c(0.18941726,0.72354824,0.08703451)
  )})
BF_fmri2b <- BF(fmri.lm, hypothesis = constraints.fmri, log=TRUE, BF.type = "AFBF")
test_that("for mlm object multiple hypotheses for same IV with log", {
  expect_equivalent(
    exp(BF_fmri2b$BFtu_confirmatory),BF_fmri2$BFtu_confirmatory
)})

#no seed required
constraints.fmri <- "Face_on_Deep = Vehicle_on_Deep < 0;
    Face_on_Middle = Vehicle_on_Middle < 0"
BF_fmri3 <- BF(x = fmri.lm, hypothesis = constraints.fmri, BF.type = "AFBF")
test_that("for mlm object multiple hypotheses on same DV", {
  expect_equivalent(
    BF_fmri3$PHP_confirmatory,c(0.05595577,0.38093401,0.56311022 )
  )})

# test multivariate t test
# via lm object
intercept <- rep(1,nrow(fmri))
lm1 <- lm(cbind(Superficial,Middle,Deep) ~ -1 + intercept, data=fmri)
BF1 <- BF(lm1,hypothesis="intercept_on_Superficial=1 & intercept_on_Middle=1.1 & intercept_on_Deep=1.36",
          BF.type = "AFBF")
BF1b <- BF(lm1,hypothesis="intercept_on_Superficial<intercept_on_Middle<intercept_on_Deep;
           intercept_on_Superficial=intercept_on_Middle=intercept_on_Deep", BF.type = "AFBF")
# via mvt_test object
mvt_fmri <- mvt_test(fmri[,3:5],null = c(1,1.1,1.36))
BF2 <- BF(mvt_fmri,hypothesis="Superficial<Middle<Deep;Superficial=Middle=Deep", BF.type = "AFBF")
test_that("test multivariate Student t", {
  expect_equivalent(
    BF1$BFtu_confirmatory,BF2$BFtu_exploratory
  )
  expect_equivalent(
    BF1b$BFtu_confirmatory,BF2$BFtu_confirmatory
  )
})

set.seed(123)
X <- as.data.frame(rmvnorm(10,mean=c(0,0)))
Y <- as.data.frame(rmvnorm(20,mean=c(1,1),sigma=.5*diag(2)+.5))
mvt_test1 <- mvt_test(X,Y,null=c(-1,-1))
BF3 <- BF(mvt_test1,hypothesis="difference_V1>difference_V2;
          difference_V1=difference_V2;
          0>difference_V1>difference_V2", BF.type = "AFBF")
test_that("test multivariate Student t, independent samples", {
  expect_equivalent(
    round(BF3$BFtu_exploratory[1],2),9.61
  )
  expect_equivalent(
    round(BF3$PHP_confirmatory[1:2],2),c(0.17,0.29)
  )
})

set.seed(123)
X <- as.data.frame(rmvnorm(10,mean=c(0,0)))
Y <- as.data.frame(rmvnorm(10,mean=c(1,1),sigma=.5*diag(2)+.5))
mvt_test1 <- mvt_test(X,Y,paired=TRUE)
BF3 <- BF(mvt_test1,hypothesis="(difference_V1,difference_V2)<0;difference_V1=difference_V2=0",log=TRUE,
          BF.type = "AFBF")
Diff1 <- X[,1] - Y[,1]
Diff2 <- X[,2] - Y[,2]
ones1 <- rep(1,length(Diff1))
mlm1 <- lm(cbind(Diff1,Diff2) ~ -1 + ones1)
BF3a <- BF(mlm1,hypothesis="(ones1_on_Diff1,ones1_on_Diff2)<0;ones1_on_Diff1=ones1_on_Diff2=0",BF.type = "AFBF")
test_that("test multivariate Student t, paired samples", {
  expect_equivalent(
    round(c(BF3$PHP_exploratory),3),c(0.519,0.481)
  )
  expect_equivalent(
    round(c(BF3$BFmatrix_confirmatory[1:2,3]),1),c(2.8,1.8)
  )
  expect_equivalent(
    round(BF3a$BFtu_confirmatory,3),round(exp(BF3$BFtu_confirmatory),3)
  )
})

set.seed(123)
X <- as.data.frame(rmvnorm(10,mean=c(0)))
mvt_test1 <- mvt_test(X)
BF3a <- BF(mvt_test1,log=TRUE, BF.type = "AFBF")
t_test1 <- t_test(X)
BF3b <- BF(t_test1,log=TRUE, BF.type = "AFBF")
test_that("test multivariate Student t, paired samples", {
  expect_equivalent(
    round(c(BF3a$BFtu_exploratory),3),round(c(BF3b$BFtu_exploratory),3)
  )
})


