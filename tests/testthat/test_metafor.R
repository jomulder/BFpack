
vi <- 1:10 #squared s.e.'s
yi <- 1:10 #sample means

### Fit a random-effects model to the data
res <- metafor::rma(yi = yi, vi = vi)
BF1 <- BF(res)

#check results exploratory test
test_that("metafor exploratory test", {
  expect_equivalent(
    round(BF1$BFtu_exploratory[,3],3),c(1.35,2), tolerance = .05
  )})
