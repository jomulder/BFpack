# test BF for t test

#####################
# ONE SAMPLE T TEST #
#####################

# test if confirmatory test gives same result exploratory test
ttest1 <- t_test(therapeutic,mu=5)
BF1 <- BF(x=ttest1,hypothesis="mu=5;mu<5",prior.hyp.conf=c(1/2,1/4,1/4),BF.type = "FBF")
test_that("1 sample t test of multiple hypotheses correctly evaluated", {
  expect_true(
    all.equal(unname(c(BF1$PHP_exploratory)),unname(BF1$PHP_confirmatory))
)})

ttest1 <- t_test(therapeutic,mu=5)
BF1 <- BF(x=ttest1,hypothesis="mu=5;mu<5",complement=FALSE, BF.type = "AFBF")
test_that("1 sample t test of multiple hypotheses correctly evaluated", {
  expect_equal(
    c(0.3523759,0.6476241),unname(BF1$PHP_confirmatory),tol=.01
  )})
BF1b <- BF(x=ttest1,hypothesis="mu=5;mu<5",complement=FALSE,prior.hyp.explo=1:3, BF.type = "AFBF")
test_that("1 sample t test of multiple hypotheses correctly evaluated", {
  expect_equal(
    unname(BF1$BFtu_exploratory * (1:3) / sum(BF1$BFtu_exploratory * (1:3))),
    unname(BF1b$PHP_exploratory),
    tol=.01
  )})
ttest1 <- t_test(therapeutic,mu=5)

BF1 <- BF(x=ttest1,hypothesis="mu=5;mu<5",complement=FALSE,log=TRUE, BF.type = "AFBF")
test_that("test log(BF) t test", {
  expect_equal(
    c(0.6086118),BF1$BFmatrix_confirmatory[2,1],tol=.01
  )})
BF1 <- BF(x=ttest1,hypothesis="mu=5;mu<6 & mu>4",complement=TRUE,log=FALSE, BF.type = "AFBF",prior.hyp.conf = c(2,1,1))
summary(BF1)
BF1$BFtable_confirmatory

# test if one-sided PMP is same as one-sided p-value
ttest2 <- t_test(therapeutic,mu=5,alternative="less")
BF2 <- BF(x=ttest2,hypothesis="mu>5", BF.type = "AFBF")
test_that("1 sample t test of one-sided hypotheses correctly evaluated", {
expect_equivalent(
  ttest2$p.value,c(BF2$PHP_confirmatory)[1]
)})

######################################
# TWO SAMPLES T TEST EQUAL VARIANCES #
######################################

# check if posterior model probabilities are correct
ttest3 <- t_test(therapeutic,therapeutic*.9+.1,var.equal=TRUE)
BF3 <- BF(ttest3,prior.hyp.explo=c(1,1,1), BF.type = "AFBF")
test_that("2 samples t test of exploratory hypotheses correctly evaluated
          with equal variances", {
  expect_equivalent(
    c(unname(BF3$PHP_exploratory)),c(0.767913,0.04941605,0.1826709)
    ,tolerance = .00001)
})
# check if posterior model probabilities are correct
ttest3a <- t_test(x=therapeutic$correct,y=therapeutic$correct*.9+.1,paired=TRUE)
BF3a <- BF(ttest3a,prior.hyp.explo=c(1,1,1), BF.type = "AFBF",hypothesis="difference=0.3;difference<0.3")
test_that("2 samples t test of exploratory hypotheses correctly evaluated
          with equal variances", {
            expect_equivalent(
              round(c(unname(BF3a$PHP_confirmatory)),2),c(0.59,0.04,0.36)
              ,tolerance = .00001)
          })

ttest3a <- t_test(x=therapeutic$correct,y=therapeutic$correct*.9+.1,paired=TRUE)
BF3a <- BF(ttest3a,prior.hyp.explo=c(1,1,1), BF.type = "FBF",hypothesis="difference < 1;difference > -1")
BF3a$BFtable_confirmatory

# t test check for testing interval hypotheses
set.seed(123)
ttest3 <- t_test(therapeutic,therapeutic*runif(length(therapeutic),min=.9,max=1.1)+.1,var.equal=TRUE)
BF3 <- BF(ttest3,hypothesis="difference< 0.5 & difference > -0.5; difference > 0.5; difference < -0.5",
          BF.type = "FBF",log = TRUE)
BF3b <- BF(ttest3,BF.type = "FBF",log = TRUE,prior.hyp.explo = 3:1)
test_that("2 samples t test of exploratory hypotheses correctly evaluated
          with equal variances", {
            expect_equivalent(
              c(unname(BF3$BFmatrix_confirmatory[1,])),c(0,3.5,4.1)
              ,tolerance = 1)
            expect_equivalent(
              BF3b$PHP_exploratory,
              exp(BF3$BFtu_exploratory)*(3:1)/sum(exp(BF3$BFtu_exploratory) * (3:1)),
              tolerance = 1)
          })

# test if one-sided PMP is same as one-sided p-value
ttest4 <- t_test(therapeutic,therapeutic*.9+.1,var.equal=TRUE,alternative="greater")
BF4 <- BF(ttest4,"difference<0", BF.type = "AFBF")
test_that("2 samples t test of one-sided hypotheses correctly evaluated
          with equal variances", {
  expect_equivalent(
    ttest4$p.value,c(BF4$PHP_confirmatory)[1]
)})


########################################
# TWO SAMPLES T TEST UNEQUAL VARIANCES #
########################################

# check posterior probabilities for a given data set
test_that("2 samples t test of two-sided hypotheses correctly evaluated
          with unequal variances", {
  skip_on_cran()
  ttest5 <- t_test(therapeutic,therapeutic*.7+2.5,"two.sided",var.equal=FALSE)
  set.seed(123)
  BF5 <- BF(ttest5,hypothesis="difference=0;difference<0",prior.hyp.explo=c(1,1,1))
  set.seed(123)
  BF5a <- BF(ttest5,prior.hyp.explo=0:2)
  expect_equivalent(
    c(unname(BF5$PHP_exploratory)),c(unname(BF5$PHP_confirmatory)),
    tolerance = .05)
  expect_equivalent(
    unname(BF5a$PHP_exploratory),
    unname(BF5$BFtu_exploratory * (0:2) / sum(BF5$BFtu_exploratory * (0:2))),
    tolerance = .05)

  BF5b <- BF(ttest5,hypothesis="difference=0; difference> -1 & difference<1; difference< -1; difference>1",
             BF.type="FBF",log=TRUE)
  expect_equivalent(
    length(unname(BF5b$PHP_confirmatory)),4
    )
  expect_equivalent(
    round(unname(BF5b$PHP_confirmatory),4),c(0.0464, 0.5373, 0.4163, 0.0000), tol=.05
  )
  BF5c <- BF(ttest5,hypothesis="difference=0; difference> -1 & difference<1; difference< -1",BF.type="FBF",log=TRUE)
  expect_equivalent(
    BF5c$BFtu_confirmatory[4],BF5b$BFtu_confirmatory[4], tol=.05
  )
})




