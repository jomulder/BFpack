# test for correlation test on hetcor object
set.seed(54)
res <- polycor::hetcor(fmri[,3:5])
#BF1 <- BF(res, hypothesis = "Deep_with_Superficial > Middle_with_Superficial")
BF1 <- BF(res)
PHPexplo <- matrix(c(0.535,  0.200,  0.265,
0.517,  0.158,  0.325,
0.538,  0.217,  0.245),nrow=3,byrow=T)
test_that("Hetcor exploratory BF gives correct result", {expect_equivalent(
  round(BF1$PHP_exploratory,3),PHPexplo)})
set.seed(463)
BF2 <- BF(res,hypothesis="(Middle_with_Superficial,Deep_with_Superficial,Deep_with_Middle) > 0;
          Middle_with_Superficial=Deep_with_Superficial=Deep_with_Middle= 0")

test_that("Hetcor two hypotheses correctly evaluated", {
  expect_equivalent(
    unname(BF2$BFtu_confirmatory)[-2],
    c(1.616808,16,0.911058)[-2],
    tolerance = .05
)
  expect_equivalent(
    unname(BF2$BFtu_confirmatory)[2],
    c(1.616808,16,0.911058)[2],
    tolerance = 1
  )})

set.seed(564)
BF3 <- BF(res,hypothesis="Middle_with_Superficial > Deep_with_Superficial")
test_that("Hetcor one order hypothesis correctly evaluated", {
  expect_equivalent(
    unname(BF3$BFtu_confirmatory),c(0.849162,1.151505), tolerance = .01
  )})

set.seed(164)
BF4 <- BF(res,hypothesis="Middle_with_Superficial = Deep_with_Superficial")
test_that("Hetcor one order hypothesis correctly evaluated", {
  expect_equivalent(
    round(unname(BF4$PHP_confirmatory),6),c(0.663753,0.336247), tolerance = .01
  )})


set.seed(164)
BF5 <- BF(res,hypothesis="Middle_with_Superficial = Deep_with_Superficial > 0")
test_that("Hetcor one hypothesis with equality and order constraint correctly evaluated", {
  expect_equivalent(
    round(unname(BF5$PHP_confirmatory),6),c(0.726939,0.273061), tolerance = .01
  )})
