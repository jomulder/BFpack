
# exploratory testing correlations in multivariate normal model
set.seed(123)
cor1 <- cor_test(mtcars[,c(1,2,3,4,6,7)], formula = ~ disp + wt)
BF1 <- BF(cor1)
PHPexplo <- matrix(
  c(0.315,  0.627,  0.059,
    0.348,  0.583,  0.069,
    0.299,  0.051,  0.651,
    0.410,  0.112,  0.478,
    0.312,  0.624,  0.064,
    0.151,  0.821,  0.027),nrow=6,byrow=T)
test_that("BF.mlm exploratory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    BF1$PHP_exploratory,PHPexplo, tolerance = .2
  )})
# confirmatory hypothesis test on the correlations
BF2 <- BF(cor1,hypothesis="cyl_with_mpg<hp_with_mpg<hp_with_cyl<0;
   cyl_with_mpg<hp_with_mpg<hp_with_cyl<0")
BFtable <- structure(c(1, 1, 1, 0.0204763291181508, 0.0208215890818035,
                       0.9815, 1, 1, 1, 0.0546186204667969, 0.0545440643456688, 0.9447,
                       1, 1, 1, 2.66740293885887, 2.6195918155611, 0.962506367804381,
                       2.3851, 2.4045, 0.962506367804381, 0.426818539062759,
                       0.419168148677558, 0.154013312259683), .Dim = c(3L, 8L), .Dimnames = list(
                         c("H1", "H2", "H3"), c("comp_E", "comp_O", "fit_E", "fit_O",
                                                "BF_E", "BF_O", "BF", "PHP")))

test_that("BF.mlm confirmatory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    as.vector(unname(BF2$BFtable_confirmatory))[-c(16, 17, 19, 20, 22, 24)], as.vector(unname(BFtable))[-c(16, 17, 19, 20, 22, 24)], tolerance = .03
  )
  expect_equivalent(
    as.vector(unname(BF2$BFtable_confirmatory))[c(22, 24)], as.vector(unname(BFtable))[c(22, 24)], tolerance = .06
  )
  expect_equivalent(
    as.vector(unname(BF2$BFtable_confirmatory))[c(16, 17, 19, 20)], as.vector(unname(BFtable))[c(16, 17, 19, 20)], tolerance = .5
  )
})








