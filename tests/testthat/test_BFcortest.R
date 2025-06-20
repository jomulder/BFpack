
# exploratory testing correlations in multivariate normal model
set.seed(123)
cor1 <- cor_test(mtcars[,4:6],iter = 1e3,burnin = 0)
BF1 <- BF(cor1,prior.hyp.explo = c(1,1,1),cov.prob=.99)
test_that("BF.cor_test use of cov.prob argument", {
  expect_equivalent(
    colnames(BF1$estimates)[3],"0.5%"
  )
  expect_equivalent(
    BF1$estimates[1,3],-0.7098442,tol=.05
  )
})
BF1a <- BF(cor1,prior.hyp.explo = 3:5)
PHPexplo <- matrix(
  c(0.14,  0.84,  0.01,
    0.0,  0.0,  1.0,
    0.0,  1.0,  0.0),nrow=3,byrow=T)
test_that("BF.cor_test exploratory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    round(BF1$PHP_exploratory,2),PHPexplo, tolerance = .1
  )
  expect_equivalent(
    BF1a$PHP_exploratory[1,],
    BF1$BFtu_exploratory[1,]*(3:5)/sum(BF1$BFtu_exploratory[1,]*(3:5)),
    tolerance = .02
  )
})
# confirmatory hypothesis test on the correlations
BF2 <- BF(cor1,hypothesis="wt_with_drat<drat_with_hp<wt_with_hp;
   wt_with_drat=drat_with_hp<0<wt_with_hp")
BFtable <- matrix(c(   0, 1.0, 5.9,
                      -1.3, 0, 4.6,
                    -5.9, -4.6, 0),byrow=TRUE,nrow=3)
test_that("BF.cor_test confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    round(log(BF2$BFmatrix_confirmatory),1),BFtable, tolerance = .5
  )})

BF2.log <- BF(cor1,hypothesis="wt_with_drat<drat_with_hp<wt_with_hp;
   wt_with_drat=drat_with_hp<0<wt_with_hp",log=TRUE)
test_that("log BF.cor_test confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    round(BF2.log$BFmatrix_confirmatory,1),BFtable, tolerance = .5
  )})

# test a single correlation
set.seed(123)
cor2 <- cor_test(mtcars[,1:2],burnin=1e2)
#print(cor2)
BF2 <- BF(cor2,hypothesis="cyl_with_mpg= -.9")
logBFexplo <- matrix(
  c(-18.46,  .69,  -22.01),nrow=1,byrow=T)
# exploratory hypothesis test on the correlation
test_that("BF.cor_test exploratory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    log(BF2$BFtu_exploratory),logBFexplo, tolerance = .3
  )})
# confirmatory hypothesis test on the correlation
test_that("BF.cor_test confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    log(BF2$BFmatrix_confirmatory[1,2]),1.99, tolerance = .2
  )})

# estimate correlations for unequal groups
set.seed(123)
cor2b <- cor_test(mtcars[1:10,2:4],mtcars[11:32,2:4])
#print(cor2b)
estimates_check <- c(.69,.64,.64,.85,.74,.67)
test_that("check estimates of correlations of multiple groups", {
  expect_equivalent(
    round(cor2b$correstimates[,1],2),estimates_check, tolerance = .1
  )})

# test a single correlation in multiple groups
test_that("BF.cor_test exploratory hypotheses on correlations correctly evaluated", {
  skip_on_cran()
  set.seed(123)
  cor3 <- cor_test(mtcars[,3:4],mtcars[,5:6])
  BF3 <- BF(cor3,hypothesis="hp_with_disp_in_g1= -wt_with_drat_in_g2")
  expect_equivalent(
    round(BF3$PHP_confirmatory,2),c(.78,.28), tolerance = .1
  )})

# test a single correlation on categorical outcomes
set.seed(123)
mtcars_test <- mtcars[,8:9]
mtcars_test[,2] <- as.factor(mtcars_test[,2])
mtcars_test[,1] <- as.factor(mtcars_test[,1])
cor2 <- cor_test(mtcars_test,burnin = 5e2,iter=3000,nugget.scale=.995)
#print(cor2)
test_that("check estimate of polychoric correlation", {
  expect_equivalent(
    round(cor2$correstimates[1,1],2),.23, tolerance = .1
  )})
BF2 <- BF(cor2,hypothesis="am_with_vs= .1",prior.hyp.explo = c(1,1,1))
PHPexplo <- matrix(
  c(.486,  .097,  .417),nrow=1,byrow=T)
# exploratory hypothesis test on the correlation
test_that("BF.cor_test exploratory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    round(BF2$PHP_exploratory,3),PHPexplo, tolerance = .1
  )})
# confirmatory hypothesis test on the correlation
test_that("BF.cor_test confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    round(log(BF2$BFmatrix_confirmatory[1,2]),2),.89, tolerance = .1
  )})

# test with combinations between continuous and ordinal (categorical) outcome variables
set.seed(187)
mtcars_test <- rbind(mtcars[,c(1,2,9,10)],mtcars[,c(1,2,9,10)])
mtcars_test[,2] <- as.ordered(mtcars_test[,2])
mtcars_test[,3] <- as.factor(mtcars_test[,3])
mtcars_test[,4] <- as.integer(mtcars_test[,4])
cor4 <- cor_test(mtcars_test,iter = 1e3,burnin = 3e2,nugget.scale = .98)
#print(cor4)
BF4 <- BF(cor4,log = TRUE)
test_that("BF.cor_test exploratory hypotheses on correlations mixed measurement levels", {
  expect_equivalent(
    round(BF4$BFtu_exploratory[,2],1),c(0.7,-13.4,-10.3,.7,.7,-18.9), tolerance = .1
  )})
#

#test ordinal correlations multiple groups
set.seed(123)
group1 <- data.frame(cbind(round(runif(20)*2+1),rnorm(20)))
class(group1$X1) <- "ordered"
group2 <- data.frame(cbind(round(runif(14)*2+1),rnorm(14)))
class(group2$X1) <- "ordered"
cor4 <- cor_test(group1,group2,iter = 1e3,burnin = 3e2)
test_that("test ordinal correlations multiple groups", {
  expect_equivalent(
    round(cor4$meanF,2),c(-0.4,0.22), tolerance = .1
  )
})

# test marginally uniform prior
# set.seed(13)
# cor5 <- cor_test(group1,group2,iter = 1e3,burnin = 3e2, prior.cor = "marg.unif")
# BF5test <- BF(cor5)
# test_that("test ordinal correlations multiple groups, marginally uniform prior", {
#   expect_equivalent(
#     round(cor4$meanF,2),c(-0.34,0.20), tolerance = .1
#   )
#   expect_equivalent(
#     round(BF5test$BFtu_exploratory[,1],3),c(1.37,2.35), tolerance = .1
#   )
# })

test_that("correlation test for intervals", {
  skip_on_cran()
  set.seed(1)
  cor_test1 <- cor_test(mtcars[,c('vs','am')])
  BF1 <- BF(x=cor_test1,hypothesis = "am_with_vs = 0; -.1<am_with_vs<.1; .1<am_with_vs<.5; -.5<am_with_vs< -.1")
  expect_equivalent(
    sum(BF1$BFtable_confirmatory[-1,2]), 1, tolerance = .05
  )
  expect_equivalent(
    BF1$BFtu_confirmatory, c(3.12, 3.16, 3.06, .38, 0.01), tolerance = .2
  )
  set.seed(1)
  BF2 <- BF(x=cor_test1,hypothesis = "am_with_vs = 0; -.1<am_with_vs<.1; -.1<am_with_vs<.5")
  expect_equivalent(
    BF2$BFtu_confirmatory, c(3.12, 3.15, 3.06, 0.13), tolerance = .2
  )
})

