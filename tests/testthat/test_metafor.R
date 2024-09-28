
vi <- 1:10 #squared s.e.'s
yi <- 1:10 #sample means

### Fit a random-effects model to the data
res <- metafor::rma(yi = yi, vi = vi)

#check results exploratory test
test_that("metafor exploratory test", {
  set.seed(123)
  BF1 <- BF(res,BF.type="correlation",iter=1e3,cov.prob=.9)
  expect_equivalent(
    round(BF1$BFtu_exploratory[,3],3),c(1.918,1.954), tolerance = .05
  )
  expect_equivalent(
    round(BF1$estimates[1:2,3],2),c(2.66,1.04), tolerance = .5
  )
})
