# test BF on lm object

lm1 <-  lm(wt ~ disp + drat + hp, mtcars)
BF1 <- BF(lm1)
#check exploratory test
PHPexplo <- matrix(
  c(0.065,  0.003,  0.931,
  0.001,  0.000,  0.999,
  0.612,  0.342,  0.046,
  0.742,  0.178,  0.080),nrow=4,byrow=T)
test_that("BF.lm exploratory hypotheses correctly evaluated", {
  expect_equivalent(
    PHPexplo,round(BF1$PHP_exploratory,3)
)})


BF2 <- BF(lm1,hypothesis="disp=drat=0;disp>drat>0;disp>drat=0")
test_that("BF.lm multiple confirmatory hypotheses correctly evaluated", {
  expect_equivalent(
    round(BF2$PHP_confirmatory,5),c(0.00000,0.53382,0.40265,0.06353)
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



