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
expect_equivalent(
  round(unname(BF2$BFtu_confirmatory),6),c(1.614833,13.796779,0.911219)
)



