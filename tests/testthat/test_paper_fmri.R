library(testthat)
constraints.fmri2 <- "Face_on_Deep = Face_on_Superficial = Face_on_Middle < 0; Face_on_Deep < Face_on_Superficial = Face_on_Middle < 0"


# Missing data ------------------------------------------------------------

fmri_missing <- fmri
set.seed(123)
for(i in 1:6){
  fmri_missing[sample(1:nrow(fmri), 1), sample(1:ncol(fmri), 1)] <- NA
}

fmri_listdel <- fmri_missing[!is.na(apply(fmri_missing, 1, sum)),]
fmri.lm2_listdel <- lm(cbind(Superficial, Middle, Deep) ~ Face + Vehicle, data = fmri_listdel)

res <- BF(fmri.lm2_listdel, hypothesis = constraints.fmri2)
test_that("BF.fmri with missing data runs", {
  expect_equivalent(res$PHP_confirmatory, c(0.12, .80, .08), tolerance = .01)
})



fmri_missing <- rbind(fmri, fmri)
set.seed(123)
for(i in 1:10){
  fmri_missing[sample(1:nrow(fmri), 1), sample(1:ncol(fmri), 1)] <- NA
}

fmri_listdel <- fmri_missing[!is.na(apply(fmri_missing, 1, sum)),]
fmri.lm2_listdel <- lm(cbind(Superficial, Middle, Deep) ~ Face + Vehicle, data = fmri_listdel)

test_that("BF.fmri with missing data throws an error", {
  expect_error(BF(fmri.lm2_listdel, hypothesis = constraints.fmri2))
})


# Without missing data ----------------------------------------------------

fmri_missing <- fmri
fmri_listdel <- fmri_missing[!is.na(apply(fmri_missing, 1, sum)),]
fmri.lm2_listdel <- lm(cbind(Superficial, Middle, Deep) ~ Face + Vehicle, data = fmri_listdel)

res <- BF(fmri.lm2_listdel, hypothesis = constraints.fmri2)
test_that("BF.fmri without missing data works", {
  expect_equivalent(res$PHP_confirmatory, c(0.05, .93, .02), tolerance = .01)
})

