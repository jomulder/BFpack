
# exploratory testing correlations in multivariate normal model
set.seed(123)
cor1 <- cor_test(mtcars[,4:6])
BF1 <- BF(cor1)
PHPexplo <- matrix(
  c(0.103,  0.888,  0.008,
    0.001,  0.000,  0.999,
    0.000,  1.000,  0.000),nrow=3,byrow=T)
test_that("BF.cor_test exploratory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    BF1$PHP_exploratory,PHPexplo, tolerance = .2
  )})
# confirmatory hypothesis test on the correlations
BF2 <- BF(cor1,hypothesis="wt_with_drat<drat_with_hp<wt_with_hp;
   wt_with_drat=drat_with_hp<0<wt_with_hp")
BFtable <- log(matrix(c(1.000, 3.495, 287.256,
                    0.286, 1.000,  82.192,
                    0.004, 0.012,   1.000),byrow=TRUE,nrow=3))
test_that("BF.cor_test confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    log(BF2$BFmatrix_confirmatory),BFtable, tolerance = .5
  )})

# confirmatory hypothesis test on the correlations
BF2 <- BF(cor1,hypothesis="wt_with_drat<drat_with_hp<wt_with_hp;
   wt_with_drat=drat_with_hp<0<wt_with_hp",complement=FALSE)
BFtable <- matrix(c(1.000, 3.446,
                    0.290, 1.000),byrow=TRUE,nrow=2)
test_that("BF.cor_test confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    unname(BF2$PHP_confirmatory),c(0.775,0.225), tolerance = .2
  )})








