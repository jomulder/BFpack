
#test for variances
vtest1 <- var_test(InsectSprays$count, InsectSprays$spray)
hypothesis <- "A=B=F>C=D=E"
set.seed(123)
BF1 <- BF(x=vtest1,hypothesis)
BF1$BFtable_confirmatory

#check results exploratory test
expect_equivalent(
  round(BF1$PHP_exploratory,7),c(0.0044175,0.9955825)
)
#check results confirmatory test
expect_equivalent(
  round(BF1$PHP_confirmatory,7),c(0.9911905,0.0088095)
)



