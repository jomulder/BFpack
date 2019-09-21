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
test_that("Hetcor two hypotheses correctly evaluated", {expect_equivalent(unname(BF2$BFtu_confirmatory), c(1.624, 14.125, 0.911), tolerance = .001)})


set.seed(564)
BF3 <- BF(res,hypothesis="Middle_with_Superficial > Deep_with_Superficial")
test_that("Hetcor two hypotheses correctly evaluated", {expect_equivalent(unname(BF2$BFtu_confirmatory), c(1.624, 14.125, 0.911), tolerance = .001)})

