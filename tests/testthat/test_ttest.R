# test BF for t test

#####################
# ONE SAMPLE T TEST #
#####################

# test if confirmatory test gives same result exploratory test
ttest1 <- bain::t_test(therapeutic,mu=5)
BF1 <- BF(x=ttest1,hypothesis="mu=5;mu<5")
test_that("1 sample t test of multiple hypotheses correctly evaluated", {
  expect_true(
    all.equal(unname(c(BF1$PHP_exploratory)),unname(BF1$PHP_confirmatory))
)})

# test if one-sided PMP is same as one-sided p-value
ttest2 <- t_test(therapeutic,mu=5,alternative="less")
BF2 <- BF(x=ttest2,hypothesis="mu>5")
test_that("1 sample t test of one-sided hypotheses correctly evaluated", {
expect_equivalent(
  ttest2$p.value,c(BF2$PHP_confirmatory)[1]
)})

######################################
# TWO SAMPLES T TEST EQUAL VARIANCES #
######################################

# check if posterior model probabilities are correct
ttest3 <- t_test(therapeutic,therapeutic*.9+.1,var.equal=T)
BF3 <- BF(ttest3)
test_that("2 samples t test of exploratory hypotheses correctly evaluated
          with equal variances", {
  expect_equivalent(
    c(unname(BF3$PHP_exploratory)),c(0.767913,0.04941605,0.1826709)
    ,tolerance = .00001)
})

# test if one-sided PMP is same as one-sided p-value
ttest4 <- t_test(therapeutic,therapeutic*.9+.1,var.equal=T,alternative="greater")
BF4 <- BF(ttest4,"difference<0")
test_that("2 samples t test of one-sided hypotheses correctly evaluated
          with equal variances", {
  expect_equivalent(
    ttest4$p.value,c(BF4$PHP_confirmatory)[1]
)})


########################################
# TWO SAMPLES T TEST UNEQUAL VARIANCES #
########################################

# check posterior probabilities for a given data set
ttest5 <- t_test(therapeutic,therapeutic*.7+2.5,"two.sided",var.equal=F)
set.seed(123)
BF5 <- BF(ttest5)
test_that("2 samples t test of two-sided hypotheses correctly evaluated
          with unequal variances", {
  expect_equivalent(
    c(unname(BF5$PHP_exploratory)),c(0.0607,0.9374,0.0019),
    tolerance = .001)
})




