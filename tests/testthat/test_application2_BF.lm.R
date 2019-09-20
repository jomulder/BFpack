aov1 <- aov(price ~ anchor*motivation,data=tvprices)
#names(aov1$coefficients)
BF1 <- BF(aov1)
#check number of tests in exploratory ANOVA
expect_equivalent(
  dim(BF1$BFtu_exploratory),c(4,3)
)
expect_equivalent(
  dim(BF1$PHP_exploratory),c(4,3)
)
expect_equivalent(
  names(aov1$coefficients),row.names(BF1$PHP_exploratory)
)
expect_equivalent(
  round(BF1$PHP_exploratory[,1],3),c(0.808,0,0,.144)
)


