cor_est <- get_estimates(cor(iris[,1:4]))

test_that("get_estimates works for correlations", {
  expect_equal(as.vector(cor_est$estimate),
               as.vector(cor(iris[,1:4])))
  })
