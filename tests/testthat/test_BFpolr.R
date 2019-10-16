options(contrasts = c("contr.treatment", "contr.poly"))
house.plr <- MASS::polr(Sat ~ Infl + Type + Cont, weights = Freq, data = MASS::housing)
house.plr

BF1 <- BF(house.plr,hypothesis="Low|Medium<0<Medium|High & TypeApartment=TypeAtrium")

test_that("BF.polr exploratory hypotheses correctly evaluated", {
  expect_equivalent(
    BF1$PHP_confirmatory,c(0.977,0.023)
  )})
