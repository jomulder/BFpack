
set.seed(57)
x <- c(rnorm(20), rnorm(20, sd = .765))
g <- factor(c(rep("V1", 20), rep("V2", 20)))
res_standard <- bartlett.test(x, g)
res <- bartlett_test(x, g)

test_that("bartlett_test same as bartlett.test", {
  expect_equivalent(res_standard, res[1:length(res_standard)])
})

#get_estimates(res)
out <- BF(res, hypothesis = "V1 > V2")

test_that("BF.bartlett returns correct results,", {
  expect_true(out$PHP_confirmatory[1] > .5 & out$PHP_confirmatory[2] < .5)
})

in_count <- c(10, 7, 20, 14, 14, 12, 10, 23, 17, 20, 14, 13, 11, 17, 21,
              11, 16, 14, 17, 17, 19, 21, 7, 13, 0, 1, 7, 2, 3, 1, 2, 1, 3,
              0, 1, 4, 3, 5, 12, 6, 4, 3, 5, 5, 5, 5, 2, 4, 3, 5, 3, 5, 3,
              6, 1, 1, 3, 2, 6, 4, 11, 9, 15, 22, 15, 16, 13, 10, 26, 26, 24,
              13)

in_spray <- structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L,
                        2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L,
                        3L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L,
                        4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 6L,
                        6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L), .Label = c("A",
           "B", "C", "D", "E", "F"), class = "factor")

res_standard <- bartlett.test(in_count, in_spray)
res <- bartlett_test(in_count, in_spray)
test_that("bartlett_test same as bartlett.test", {
  expect_equivalent(res_standard, res[1:length(res_standard)])
})

out <- BF(res, "F> A=B>D>E=C")
test_that("BF.bartlett returns correct results,", {
  expect_equivalent(out$PHP_confirmatory, c(.9958, .00417), tolerance = .001)
})

set.seed(57)
out2 <- BF(res, "F> A=B>D>E=C", log = TRUE)
test_that("BF.bartlett returns correct results,", {
  expect_equivalent(out2$BFmatrix_confirmatory[2,1], -5.452534, tolerance = .01)
})

set.seed(57)
out2b <- BF(res, log = TRUE, prior.hyp.explo = c(1,2))
test_that("BF.bartlett returns correct results with prior.hyp.explo,", {
  expect_equivalent(out2b$PHP_exploratory,
                    exp(out2$BFtu_exploratory) * c(1,2) / sum(exp(out2$BFtu_exploratory) * c(1,2)),
                    tolerance = .01)
})




