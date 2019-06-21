# test BF for t test

# one sample t test
library(bain)
ttest1 <- t_test(therapeutic,mu=5)
BF(ttest1)
BF(ttest1,"mu=5;mu<5")

# two samples t test with equal variances
ttest1 <- t_test(therapeutic,therapeutic*.7+.5,var.equal=T)
BF(ttest1)
BF(ttest1,hypothesis="difference=0.2;difference<0.2")

# two samples t test with unequal variances
ttest1 <- t_test(therapeutic,therapeutic*.7+.5,"two.sided",var.equal=F)
BF(ttest1)
BF(ttest1,hypothesis="difference=0.2;difference<0.2")





