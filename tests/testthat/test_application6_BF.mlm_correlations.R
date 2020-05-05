set.seed(123)
#split data frame in two data frame for each group
memoryHC <- subset(memory,Group=="HC")[,-7]
memorySZ <- subset(memory,Group=="SZ")[,-7]
cor1 <- cor_test(memoryHC,memorySZ)
BF6_cor <- BF(cor1, hypothesis=
    "Del_with_Im_in_1 > Del_with_Im_in_2 &
    Del_with_Wmn_in_1 > Del_with_Wmn_in_2 &
    Del_with_Cat_in_1 > Del_with_Cat_in_2 &
    Del_with_Fas_in_1 > Del_with_Fas_in_2 &
    Del_with_Rat_in_1 > Del_with_Rat_in_2 &
    Im_with_Wmn_in_1 > Im_with_Wmn_in_2 &
    Im_with_Cat_in_1 > Im_with_Cat_in_2 &
    Im_with_Fas_in_1 > Im_with_Fas_in_2 &
    Im_with_Rat_in_1 > Im_with_Rat_in_2 &
    Wmn_with_Cat_in_1 > Wmn_with_Cat_in_2 &
    Wmn_with_Fas_in_1 > Wmn_with_Fas_in_2 &
    Wmn_with_Rat_in_1 > Wmn_with_Rat_in_2 &
    Cat_with_Fas_in_1 > Cat_with_Fas_in_2 &
    Cat_with_Rat_in_1 > Cat_with_Rat_in_2 &
    Fas_with_Rat_in_1 > Fas_with_Rat_in_2")
#check results
test_that("correlation test on cor_test object with two groups correctly evaluated", {
  expect_equivalent(
    log(BF6_cor$BFmatrix_confirmatory[1,2]), 8.7, tolerance = 1
)})


