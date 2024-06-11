fit <- glm(sent ~ ztrust + zfWHR + zAfro + glasses + attract + maturity +
             tattoos, family = binomial(), data = wilson)
BF_glm <- BF(fit)

ct <- lmtest::coeftest(fit)
BFdefault1 <- BF(ct[,1], Sigma = diag(ct[,2]^2), n = attr(ct,"nobs"))
BFcoeftest1 <- BF(ct)

test_that("BF.coeftest gives same exploratory results as BF.default", {
  expect_equivalent(
    BFdefault1$BFtu_exploratory,BFcoeftest1$BFtu_exploratory
  )})


