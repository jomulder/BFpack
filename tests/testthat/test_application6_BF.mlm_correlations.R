set.seed(123)
#split data frame in two data frame for each group
memoryHC <- subset(memory,Group=="HC")[,-(4:7)]
memorySZ <- subset(memory,Group=="SZ")[,-(4:7)]
cor1 <- cor_test(memoryHC,memorySZ)
BF6_cor <- BF(cor1, hypothesis=
    "Del_with_Im_in_g1 > Del_with_Im_in_g2 &
    Del_with_Wmn_in_g1 > Del_with_Wmn_in_g2 &
    Im_with_Wmn_in_g1 > Im_with_Wmn_in_g2")
#check results
test_that("correlation test on cor_test object with two groups correctly evaluated", {
  expect_equivalent(
    log(BF6_cor$BFmatrix_confirmatory[1,2]), 5.0, tolerance = .5
)})







