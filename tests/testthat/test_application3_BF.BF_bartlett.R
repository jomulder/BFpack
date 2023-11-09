bartlett <- bartlett_test(x = attention$accuracy, g = attention$group)
hypothesis <- c("Controls=TS<ADHD;
    Controls<TS=ADHD;
    Controls=TS=ADHD")
set.seed(358)
BF_var <- BF(bartlett, hypothesis)
#check posterior probability exploratory test
test_that("bartlett_test works for an exploratory test of equality of variances", {
expect_equivalent(
  round(BF_var$PHP_exploratory[1],7),0.8032805
)})

#check equality part equals NA
test_that("bartlett_test works for an exploratory test of equality of variances", {
expect_true(
  is.na(BF_var$BFtable_confirmatory[1,1])
)
})

#check names hypotheses confimatory test
test_that("names match of output for a bartlett_test", {
expect_equivalent(
  colnames(BF_var$BFmatrix_confirmatory),BF_var$hypotheses
)})
test_that("names match of output for a bartlett_test", {
expect_equivalent(
  row.names(BF_var$BFmatrix_confirmatory),BF_var$hypotheses
)})
test_that("bartlett_test works for a specific conformator test of equality of variances", {
expect_equivalent(
  round(BF_var$PHP_confirmatory,3),c(0.426,0.278,0.238,0.058)
)})

BF_var2 <- BF(bartlett, hypothesis, log = TRUE)
test_that("names match of output for a bartlett_test", {
  expect_equivalent(
    round(BF_var$BFmatrix_confirmatory,2),round(exp(BF_var2$BFmatrix_confirmatory),2),tol=.02
  )})



