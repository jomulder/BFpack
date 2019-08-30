# test BF on lm object

lm1 <-  lm(wt ~ disp + drat + hp, mtcars)
BF1 <- BF(lm1)
#check exploratory test
PHPexplo <- matrix(
  c(0.065,  0.003,  0.931,
  0.001,  0.000,  0.999,
  0.612,  0.342,  0.046,
  0.742,  0.178,  0.080),nrow=4,byrow=T)
expect_equivalent(
  PHPexplo,round(BF1$PHP_exploratory,3)
)
#check confirmatory tes
BF2 <- BF(lm1,hypothesis="disp=drat=0;disp>drat>0;disp>drat=0")
expect_equivalent(
  round(BF2$PHP_confirmatory,5),c(0.00000,0.53382,0.40265,0.06353)
)

