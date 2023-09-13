# test for correlation test on hetcor object
set.seed(54)
res <- polycor::hetcor(fmri[,3:5])
#BF1 <- BF(res, hypothesis = "Deep_with_Superficial > Middle_with_Superficial")

BF2 <- BF(res,hypothesis="(Superficial_with_Middle,Deep_with_Superficial,Deep_with_Middle) > 0;
          Middle_with_Superficial=Deep_with_Superficial=Deep_with_Middle= 0")
PHPexplo <- matrix(c(0.51,  0.21,  0.27,
0.50,  0.17,  0.33,
0.51,  0.22,  0.26),nrow=3,byrow=T)
test_that("Hetcor exploratory BF gives correct result", {expect_equivalent(
  BF2$PHP_exploratory, PHPexplo, tolerance = .01)})
set.seed(463)

test_that("Hetcor two hypotheses correctly evaluated", {
  expect_equivalent(
    unname(BF2$BFtu_confirmatory)[-2],
    c(1.63,0.91),
    tolerance = .05
)
  expect_equivalent(
    unname(BF2$BFtu_confirmatory)[2],
    c(1.616808,16,0.911058)[2],
    tolerance = 1
  )}
)

set.seed(164)
BF5 <- BF(res,hypothesis="Middle_with_Superficial = Deep_with_Superficial > 0")
test_that("Hetcor one hypothesis with equality and order constraint correctly evaluated", {
  expect_equivalent(
    unname(BF5$PHP_confirmatory),c(0.75,0.24), tolerance = .1
  )})


set.seed(564)
res1 <- polycor::hetcor(fmri[,3:4])
BF3 <- BF(res1,hypothesis="Superficial_with_Middle > .4;Superficial_with_Middle = .4; Superficial_with_Middle < .4")
test_that("Hetcor test BF3 (automatically omit complement)", {
  expect_equivalent(
    log(unname(BF3$BFtu_confirmatory)),c(-1.2,.38,.27), tolerance = .5
  )})


