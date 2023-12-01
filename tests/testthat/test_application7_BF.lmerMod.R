
#check results confirmatory test
test_that("lmerMod two hypotheses correctly evaluated", {
  skip_on_cran()
  library(lme4)
  timssICC_subset <- timssICC[(timssICC$groupNL11==1)+(timssICC$groupHR11==1)>0,
                              ][c(1:150,1395+1:150),]
  outlme1 <- lme4::lmer(math ~ -1 + lln +
                          groupNL11 + (0+groupNL11 | schoolID) +
                          groupHR11 + (0+groupHR11 | schoolID),
                        data=timssICC_subset)
  set.seed(123)
  BFicc <- BF(outlme1,hypothesis=
                "groupNL11<groupHR11;
               groupNL11=groupHR11",log=TRUE)
  expect_equivalent(
    round(BFicc$PHP_confirmatory,7),c(0.7,0.25,0.03), tolerance = .1
  )
  expect_equivalent(
    round(BFicc$BFmatrix_confirmatory[,1],7),c(0,-.93,-3), tolerance = .2
  )
  set.seed(123)
  BFicc2 <- BF(outlme1,log=TRUE,prior.hyp.explo = 3:5)
  expect_equivalent(
    round(BFicc2$PHP_exploratory[1,],7),
    round(exp(BFicc2$BFtu_exploratory[1,])*(3:5)/sum(exp(BFicc2$BFtu_exploratory[1,])*(3:5)),7),
    tolerance = .2
  )
})



#this is a smaller subset of the data with unbalanced groups
#check results exploratory test
test_that("lmerMod exploratory unbalanced", {
  skip_on_cran()
  library(lme4)
  timssICC_subset <- timssICC[timssICC$groupNL11==1,][c(1:15,15+1:12,30+1:8,45+1:4,60+1:9),]
  outlme2 <- lme4::lmer(math ~ -1 + gender +
                          (0+groupNL11 | schoolID),
                        data=timssICC_subset)
  set.seed(123)
  BFicc2 <- BF(outlme2)
  expect_equivalent(
    round(BFicc2$BFtu_exploratory,3),c(0,0,1.07), tolerance = .1
  )
})



