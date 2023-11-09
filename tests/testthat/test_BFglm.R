#Test glm class object
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
glm.D93 <- glm(counts ~ treatment, family = poisson())
set.seed(432)
BF1 <- BF(x=glm.D93, hypothesis = "treatment2 = 0; treatment2 < 0")
#check if confirmatory and exploratory test are the same for Gaussian estimator
test_that("BF.glm two hypotheses correctly evaluated", {
expect_true(
  all.equal(unname(BF1$PHP_exploratory[2,]),
          unname(BF1$PHP_confirmatory), tolerance = .005)
)})

BF1.log <- BF(x=glm.D93, hypothesis = "treatment2 = 0; treatment2 < 0", log = TRUE)
#check if confirmatory and exploratory test are the same for Gaussian estimator
test_that("BF.glm two hypotheses correctly evaluated", {
  expect_true(
    all.equal(BF1$BFmatrix_confirmatory,
              exp(BF1.log$BFmatrix_confirmatory), tolerance = .005)
  )})

