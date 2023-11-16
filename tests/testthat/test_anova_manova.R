
set.seed(123)
n <- 40
data.example <- data.frame(group1 = round(runif(n)*2), group2 = round(runif(n)), dv1 = rnorm(n), dv2 = rnorm(n,sd=.5),cov=round(rnorm(n,sd=.5)))
data.example$group1 <- as.factor(data.example$group1)
data.example$group2 <- as.factor(data.example$group2)

m1 <- aov(dv1 ~ 1 + group1, data.example)
BF1 <- BF(m1)
m2 <- aov(dv1 ~ -1 + group1, data.example)
BF2 <- BF(m2)
test_that("BF.aov exploratory ANOVA (1 factor) when including/excluding the intercept", {
  expect_equivalent(
    BF1$BFtu_exploratory,BF2$BFtu_exploratory
  )})

m3 <- aov(dv1 ~ 1 + group1 + group2, data.example)
BF3 <- BF(m3)
m4 <- aov(dv1 ~ -1 + group1 + group2, data.example)
BF4 <- BF(m4)
test_that("BF.aov exploratory ANOVA (2 factors) when including/excluding the intercept", {
  expect_equivalent(
    BF3$BFtu_exploratory,BF4$BFtu_exploratory
  )})

test_that("BF.manova exploratory MANOVA (1 factor)", {
  skip_on_cran()
  m5 <- manova(cbind(dv1,dv2) ~ 1 + group1, data.example)
  m6 <- manova(cbind(dv1,dv2) ~ -1 + group1, data.example)
  set.seed(123)
  BF5 <- BF(m5)
  set.seed(123)
  BF6 <- BF(m6)
  expect_equivalent(
    BF5$BFtu_exploratory,BF6$BFtu_exploratory
  )
})

test_that("BF.manova exploratory MANOVA (2 factors)", {
  skip_on_cran()
  m5 <- manova(cbind(dv1,dv2) ~ 1 + group1 + group2, data.example)
  m6 <- manova(cbind(dv1,dv2) ~ -1 + group1 + group2, data.example)
  set.seed(123)
  BF5 <- BF(m5)
  set.seed(123)
  BF6 <- BF(m6)
  expect_equivalent(
    BF5$BFtu_exploratory,BF6$BFtu_exploratory
  )
})

m7 <- aov(dv1 ~ 1 + group1 * group2, data.example)
BF7 <- BF(m7)
m8 <- aov(dv1 ~ -1 + group1 * group2, data.example)
BF8 <- BF(m8)
test_that("BF.aov exploratory ANOVA (2 factors w/ interaction) when including/excluding the intercept", {
  expect_equivalent(
    BF7$BFtu_exploratory,BF8$BFtu_exploratory
  )})

test_that("BF.manova exploratory MANOVA (2 factors w/ interaction)", {
  skip_on_cran()
  m5 <- manova(cbind(dv1,dv2) ~ 1 + group1 * group2, data.example)
  m6 <- manova(cbind(dv1,dv2) ~ -1 + group1 * group2, data.example)
  set.seed(123)
  BF5 <- BF(m5)
  set.seed(123)
  BF6 <- BF(m6)
  expect_equivalent(
    BF5$BFtu_exploratory,BF6$BFtu_exploratory
  )
})

test_that("BF.manova exploratory MANOVA (2 factors w/ interaction and numeric covariate, group identifier)", {
  skip_on_cran()
  m5 <- manova(cbind(dv1,dv2) ~ 1 + group1 * group2, data.example)
  m6 <- manova(cbind(dv1,dv2) ~ -1 + group1 * group2 + cov, data.example)
  set.seed(123)
  BF5 <- BF(m5)
  set.seed(123)
  BF6 <- BF(m6)
  expect_equivalent(
    BF5$fraction_group_identifier,BF6$fraction_group_identifier
  )
})

test if interaction effects are properly tested (check how RE looks like)


