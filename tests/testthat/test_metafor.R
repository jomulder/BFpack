
vi <- 1:10 #squared s.e.'s
yi <- 1:10 #sample means

### Fit a random-effects model to the data
res <- metafor::rma(yi = yi, vi = vi)

#check results exploratory test
test_that("metafor exploratory test", {
  skip_on_cran()
  set.seed(123)
  BF1 <- BF(res,BF.type="correlation")
  expect_equivalent(
    round(BF1$BFtu_exploratory[,3],3),c(1.910,1.924), tolerance = .05
  )
  expect_equivalent(
    round(BF1$estimates[1:2,1],2),c(2.14,15.41), tolerance = .5
  )
})
