# test for correlation test on hetcor object
res <- polycor::hetcor(fmri[,3:5])
BF1 <- BF(res)
PHPexplo <- matrix(c(0.535,  0.200,  0.265,
0.517,  0.158,  0.325,
0.538,  0.217,  0.245),nrow=3,byrow=T)
expect_equivalent(
  round(BF1$PHP_exploratory,3),PHPexplo
)
set.seed(123)
BF2 <- BF(res,hypothesis="(Middle_with_Superficial,Deep_with_Superficial,Deep_with_Middle) > 0;
          Middle_with_Superficial=Deep_with_Superficial=Deep_with_Middle= 0")
test_that("Hetcor two hypotheses correctly evaluated", {
  expect_equivalent(
    round(unname(BF2$BFtu_confirmatory),6),c(1.621700,14.558565,0.910662)
)})

set.seed(564)
BF3 <- BF(res,hypothesis="Middle_with_Superficial > Deep_with_Superficial")
test_that("Hetcor one order hypothesis correctly evaluated", {
  expect_equivalent(
    round(unname(BF3$BFtu_confirmatory),6),c(0.849162,1.151505)
  )})

set.seed(164)
BF4 <- BF(res,hypothesis="Middle_with_Superficial = Deep_with_Superficial")
test_that("Hetcor one order hypothesis correctly evaluated", {
  expect_equivalent(
    round(unname(BF4$PHP_confirmatory),6),c(0.663753,0.336247)
  )})

set.seed(164)
BF5 <- BF(res,hypothesis="Middle_with_Superficial = Deep_with_Superficial > 0")
test_that("Hetcor one hypothesis with equality and order constraint correctly evaluated", {
  expect_equivalent(
    round(unname(BF5$PHP_confirmatory),6),c(0.726939,0.273061)
  )})


