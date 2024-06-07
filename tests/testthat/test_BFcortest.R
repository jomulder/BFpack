
# exploratory testing correlations in multivariate normal model
set.seed(123)
cor1 <- cor_test(mtcars[,4:6],iter = 3e3,burnin = 1e3)
BF1 <- BF(cor1)
BF1a <- BF(cor1,prior.hyp.explo = 3:5)
PHPexplo <- matrix(
  c(0.08,  0.92,  0.0,
    0.0,  0.0,  1.0,
    0.0,  1.0,  0.0),nrow=3,byrow=T)
test_that("BF.cor_test exploratory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    BF1$PHP_exploratory,PHPexplo, tolerance = .1
  )
  expect_equivalent(
    BF1a$PHP_exploratory[1,],
    BF1$BFtu_exploratory[1,]*(3:5)/sum(BF1$BFtu_exploratory[1,]*(3:5)),
    tolerance = .02
  )
})
# confirmatory hypothesis test on the  correlations
BF2 <- BF(cor1,hypothesis="wt_with_drat<drat_with_hp<wt_with_hp;
   wt_with_drat=drat_with_hp<0<wt_with_hp")
BFtable <- matrix(c(   0, 1.1, 5.62,
                      -1.1, 0, 4.6,
                    -5.62, -4.6, 0),byrow=TRUE,nrow=3)
test_that("BF.cor_test confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    log(BF2$BFmatrix_confirmatory),BFtable, tolerance = .5
  )})

BF2.log <- BF(cor1,hypothesis="wt_with_drat<drat_with_hp<wt_with_hp;
   wt_with_drat=drat_with_hp<0<wt_with_hp",log=TRUE)
test_that("log BF.cor_test confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    BF2.log$BFmatrix_confirmatory,BFtable, tolerance = .5
  )})

# test a single correlation
set.seed(123)
cor2 <- cor_test(mtcars[,1:2],burnin=1e2)
BF2 <- BF(cor2,hypothesis="cyl_with_mpg= -.9")
logBFexplo <- matrix(
  c(-21.45,  .69,  -25.07),nrow=1,byrow=T)
# exploratory hypothesis test on the correlation
test_that("BF.cor_test exploratory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    log(BF2$BFtu_exploratory),logBFexplo, tolerance = .3
  )})
# confirmatory hypothesis test on the correlation
test_that("BF.cor_test confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    log(BF2$BFmatrix_confirmatory[1,2]),2.5, tolerance = .2
  )})

# estimate correlations for unequal groups
set.seed(123)
cor2b <- cor_test(mtcars[1:10,2:4],mtcars[11:32,2:4])
estimates_check <- c(.85,.81,.8,.9,.81,.76)
test_that("check estimates of correlations of multiple groups", {
  expect_equivalent(
    cor2b$correstimates[,1],estimates_check, tolerance = .1
  )})

# test a single correlation in multiple groups
test_that("BF.cor_test exploratory hypotheses on correlations correctly evaluated", {
  skip_on_cran()
  set.seed(123)
  cor3 <- cor_test(mtcars[,3:4],mtcars[,5:6])
  BF3 <- BF(cor3,hypothesis="hp_with_disp_in_g1= -wt_with_drat_in_g2")
  expect_equivalent(
    BF3$PHP_confirmatory,c(.77,.23), tolerance = .1
  )})

# test a single correlation on categorical outcomes
set.seed(123)
mtcars_test <- mtcars[,8:9]
mtcars_test[,2] <- as.factor(mtcars_test[,2])
mtcars_test[,1] <- as.factor(mtcars_test[,1])
cor2 <- cor_test(mtcars_test,burnin = 0, iter = 20)
print(cor2$meanF)
print(cor2$covmF)
print(cor2$correstimates)
print(head(cor2$corrdraws[[1]][,2,1]))
print(tail(cor2$corrdraws[[1]][,2,1]))
print(head(cor2$res$BDrawsStore[,1,1,]))
print(tail(cor2$res$BDrawsStore[,1,1,]))
print(head(cor2$res$sigmaDrawsStore[,1,]))
print(tail(cor2$res$sigmaDrawsStore[,1,]))
print(head(cor2$res$gLiuSab[,1,]))
print(tail(cor2$res$gLiuSab[,1,]))
print("W2")
print(head(cor2$res$WgroupsStore[2,1,,]))
print("W3")
print(head(cor2$res$WgroupsStore[3,1,,]))
print(head(cor2$res$CheckStore[1,1,,1,],2))
print(head(cor2$res$CheckStore[1,1,,2,],2))
print(head(cor2$res$CheckStore[2,1,,1,],2))
print(head(cor2$res$CheckStore[2,1,,2,],2))
print(head(cor2$res$CheckStore[3,1,,1,],2))
print(head(cor2$res$CheckStore[3,1,,2,],2))

test_that("check estimate of polychoric correlation", {
  expect_equivalent(
    cor2$correstimates[1,1],.2, tolerance = .1
  )})
BF2 <- BF(cor2,hypothesis="am_with_vs= .1")
PHPexplo <- matrix(
  c(.51,  .046,  .45),nrow=1,byrow=T)
# exploratory hypothesis test on the correlation
test_that("BF.cor_test exploratory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    BF2$PHP_exploratory,PHPexplo, tolerance = .1
  )})
# confirmatory hypothesis test on the correlation
test_that("BF.cor_test confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    log(BF2$BFmatrix_confirmatory[1,2]),1.35, tolerance = .1
  )})

# test with combinations between continuous and ordinal (categorical) outcome variables
set.seed(187)
mtcars_test <- rbind(mtcars[,c(1,2,9,10)],mtcars[,c(1,2,9,10)])
mtcars_test[,2] <- as.ordered(mtcars_test[,2])
mtcars_test[,3] <- as.factor(mtcars_test[,3])
mtcars_test[,4] <- as.integer(mtcars_test[,4])
cor4 <- cor_test(mtcars_test,iter = 3e3,burnin = 1e3)
BF4 <- BF(cor4,log = TRUE)
test_that("BF.cor_test exploratory hypotheses on correlations mixed measurement levels", {
  expect_equivalent(
    BF4$BFtu_exploratory[,2],c(0.7,-7.4,-9.2,.7,.7,-7), tolerance = .1
  )})

