# simulated data with correlated covariates that have strong effects
set.seed(1)
n <- 200
x1 <- rnorm(n)
x2 <- 0.7*x1 + rnorm(n, sd = .7)
x3 <- rbinom(n, 1, .5) + 1
y1 <- 1*x1 + .5*x2 + .5*x3 + rnorm(n)
y2 <- .3*y1 + .8*x1 + rnorm(n)
dat <- data.frame(y1, y2, x1, x2, x3)

test_that("cor_test results do not depend on the order of the covariates", {
  set.seed(123)
  fit1 <- cor_test(dat, formula = ~ x1 + x2 + x3, iter = 3000, burnin = 500)
  set.seed(123)
  fit2 <- cor_test(dat, formula = ~ x3 + x2 + x1, iter = 3000, burnin = 500)
  expect_equal(fit1$correstimates[, "median"], fit2$correstimates[, "median"], tolerance = .02)
})

test_that("posterior sds of the regression coefficients match lm standard errors", {
  set.seed(123)
  fit <- cor_test(dat, formula = ~ x1 + x2 + x3, iter = 3000, burnin = 500)
  B_sd <- apply(fit$res[[1]]$BDrawsStore[, 1, , ], c(2, 3), sd)
  dat_s <- as.data.frame(scale(dat))
  fit_lm <- lm(cbind(y1, y2) ~ x1 + x2 + x3, data = dat_s)
  se_lm <- matrix(sqrt(diag(vcov(fit_lm))), ncol = 2)
  expect_equal(c(B_sd), c(se_lm), tolerance = .15)
})

test_that("cor_test works when a covariate column precedes the dependent variables", {
  set.seed(123)
  expect_no_error(fit <- cor_test(mtcars[, 1:3], formula = ~ mpg, iter = 1000, burnin = 100))
  expect_equal(rownames(fit$correstimates), "disp_with_cyl")
})
