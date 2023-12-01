# test BF on lm object

lm1 <-  lm(wt ~ disp + drat + hp, mtcars)
BF1a <- BF(lm1,BF.type = 2)
BF1aa <- BF(lm1,BF.type = 2,prior.hyp.explo = 1:3)
#check exploratory test
PHPexplo <- matrix(
  c(0.065,  0.003,  0.931,
  0.001,  0.000,  0.999,
  0.612,  0.342,  0.046,
  0.742,  0.178,  0.080),nrow=4,byrow=T)
test_that("BF.lm exploratory hypotheses correctly evaluated", {
  expect_equivalent(
    PHPexplo,round(BF1a$PHP_exploratory,3)
  )
  expect_equivalent(
    unname(round(BF1aa$PHP_exploratory[3,],3)),
           unname(round(BF1a$BFtu_exploratory[3,] * c(1:3) / sum(BF1a$BFtu_exploratory[3,] * c(1:3)),3))
  )
  expect_equivalent(
    unname(round(BF1aa$PHP_exploratory[4,],3)),
    unname(round(BF1a$BFtu_exploratory[4,] * c(1:3) / sum(BF1a$BFtu_exploratory[4,] * c(1:3)),3))
  )
})

BF1b <- BF(lm1,BF.type = 1)
PHPexplo <- matrix(
  c(0.107,  0.006,  0.887,
    0.002,  0.000,  0.998,
    0.646,  0.300,  0.054,
    0.748,  0.168,  0.085),nrow=4,byrow=T)
test_that("BF.lm exploratory hypotheses correctly evaluated", {
  expect_equivalent(
    PHPexplo,round(BF1b$PHP_exploratory,3)
  )})


BF2 <- BF(lm1,hypothesis="disp=drat=0;disp>drat>0;disp>drat=0")
test_that("BF.lm multiple confirmatory hypotheses correctly evaluated", {
  expect_equivalent(
    round(BF2$PHP_confirmatory,5),c(0.000,0.533,0.402,0.064), tolerance = .01
)})

BF2 <- BF(lm1,hypothesis="disp=drat=0;disp>drat>0;disp>drat=0",complement=FALSE)
test_that("BF.lm multiple confirmatory hypotheses correctly evaluated", {
  expect_equivalent(
    round(BF2$PHP_confirmatory,5),c(0.00000,0.57003,0.42997)
  )})

BF3 <- BF(lm1,hypothesis="disp=hp=0")
test_that("BF.lm one equality hypothesis correctly evaluated", {
  expect_equivalent(
    round(BF3$PHP_confirmatory,5),c(0.00003,0.99997)
  )})

BF4 <- BF(lm1,hypothesis="drat<hp<disp")
test_that("BF.lm one order hypothesis correctly evaluated", {
  expect_equivalent(
    round(BF4$PHP_confirmatory,5),c(0.96872,0.03128)
  )})

BF5 <- BF(lm1,hypothesis="drat<hp=disp")
test_that("BF.lm one equal/order hypothesis correctly evaluated", {
  expect_equivalent(
    round(BF5$PHP_confirmatory,5),c(0.47961,0.52039)
  )})

BF5b <- BF(lm1,hypothesis="drat<hp=disp",BF.type = 1)
test_that("BF.lm one equal/order hypothesis correctly evaluated", {
  expect_equivalent(
    round(BF5b$PHP_confirmatory,5),c(0.44049,0.55951)
  )})



