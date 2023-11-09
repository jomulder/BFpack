
#test for variances
vtest1 <- bartlett_test(InsectSprays$count, InsectSprays$spray)
hypothesis <- "A=B=F>C=D=E"
set.seed(123)
BF1 <- BF(x=vtest1,hypothesis)

#check results exploratory test
test_that("BF.bartlett_htest exploratory hypotheses correctly evaluated", {
  expect_equivalent(
    round(BF1$PHP_exploratory,7),c(0.0044175,0.9955825)
)})
#check results confirmatory test
test_that("BF.bartlett_htest confirmatory hypotheses correctly evaluated", {
expect_equivalent(
  round(BF1$PHP_confirmatory,7),c(0.9911905,0.0088095)
)})

hypothesis <- "A=B=F>C=D=E; A=B=F>C>D>E"
set.seed(123)
BF1 <- BF(x=vtest1,hypothesis,complement = F, log = TRUE)
#check results confirmatory test
test_that("BF.bartlett_htest confirmatory hypotheses correctly evaluated log(BF)", {
  expect_equivalent(
    round(BF1$BFtu_confirmatory,3),c(4.723,4.245)
  )})
test_that("BF.bartlett_htest confirmatory hypotheses correctly evaluated log(BF)", {
  expect_equivalent(
    round(BF1$BFtu_exploratory[1],4),-5.4178
  )})

