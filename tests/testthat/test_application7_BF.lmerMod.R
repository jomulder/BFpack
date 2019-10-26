#this is a smaller subset of the data to decrease the computation time
timssICC_subset <- timssICC[(timssICC$groupNL11==1)+(timssICC$groupHR11==1)>0,][c(1:150,1395+1:150),]
outlme1 <- lme4::lmer(math ~ -1 + gender + weight + lln +
                  groupNL11 + (0+groupNL11 | schoolID) +
                  groupHR11 + (0+groupHR11 | schoolID),
                data=timssICC_subset)
set.seed(123)
BFicc <- BF(outlme1,hypothesis=
              "groupNL11<groupHR11;
               groupNL11=groupHR11")

#check results confirmatory test
test_that("lmerMod two hypotheses correctly evaluated", {
  expect_equivalent(
    round(BFicc$PHP_confirmatory,7),c(0.593,0.234,0.173), tolerance = .01
)})
#check results exploratory test
matrixPHP <- matrix(c(0.845, 0.052, 0.103,
                      0.000, 0.000, 1.000),nrow=2,byrow=T)
test_that("lmerMod explortatory hypotheses correctly evaluated", {
  expect_equivalent(
    BFicc$PHP_exploratory,matrixPHP
)})


