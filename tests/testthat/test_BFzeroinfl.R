fm_zip <- pscl::zeroinfl(art ~ . | 1, data = pscl::bioChemists)
set.seed(752)
BF1 <- BF(fm_zip)
test_that("BF.zeroinfl exploratory hypotheses correctly evaluated", {
  expect_equivalent(
    BF1$PHP_exploratory[,1],c(0.000,0.006,0.674,0.006,0.938,0.000,0.000)
)})

