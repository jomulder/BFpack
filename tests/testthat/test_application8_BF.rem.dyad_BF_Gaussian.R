
CovEventEff <- array(NA, dim = c(3, nrow(actors), nrow(actors)))
CovEventEff[1,,] <- as.matrix(same_division)
CovEventEff[2,,] <- as.matrix(same_hierarchy)
CovEventEff[3,,] <- as.matrix(same_building)
set.seed(9227)
fit <- relevent::rem.dyad(edgelist = relevents, n = nrow(actors),
                          effects = "CovEvent", ordinal = FALSE,
                          covar = list(CovEvent = CovEventEff), hessian = TRUE,
                          fit.method = "BPM")
# Set effect names
names(fit$coef) <- c("division", "hierarchy", "building")
# Define the hypotheses
hyp <- "division = hierarchy = building;
    division > hierarchy = building;
    division > hierarchy > building;
    division > building > hierarchy"
# Run BF
set.seed(8389)
BF_rem <- BF(x = fit, hypothesis = hyp)
# check results
test_that("BF.rem.dyad four hypotheses correctly evaluated", {
  expect_equivalent(
    BF_rem$PHP_confirmatory,c(0.000,0.715,0.063,0.222,0.000)
)})

