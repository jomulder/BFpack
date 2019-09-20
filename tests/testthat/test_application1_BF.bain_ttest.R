
# Frequentist one sample t test at the null point mu=5
ttest1 <- t_test(therapeutic,alternative="greater",mu=5)
# one-sided Bayesian one sample t test
BF1 <- BF(ttest1,hypothesis="mu<5",prior=c(.5,.5))
expect_equivalent(
  round(ttest1$p.value,3),round(BF1$PHP_confirmatory[1],3)
)
# confirmatory Bayesian one sample t test
BF2 <- BF(ttest1,hypothesis="mu=5;mu>5",prior=c(.5,.5,0))
#check prior output
expect_equivalent(
  BF2$prior,c(.5,.5,0)
)
#check posterior probabilities
expect_equivalent(
  round(BF2$PHP_confirmatory,3),c(0.943,0.057,0.000)
)



