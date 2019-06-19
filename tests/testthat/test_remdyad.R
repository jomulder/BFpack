# Load packages
library(relevent)				# Fitting a relational event model
library(bain)						# Bayes factors for informative hypotheses

# Step 1: Fit a relational event model
# ---------------------------------------
# Load package
library(relevent)

# Combine the event statistics into a CovEvent array:
# a 3-dimensional p by n by n array
CovEventEff <- array(NA, dim = c(3, nrow(actors), nrow(actors)))
CovEventEff[1,,] <- as.matrix(same_division)
CovEventEff[2,,] <- as.matrix(same_hierarchy)
CovEventEff[3,,] <- as.matrix(same_building)

# Fit a relational event model
set.seed(9227)
fit <- rem.dyad(edgelist = relevents, n = nrow(actors), effects = "CovEvent",
                ordinal = FALSE, covar = list(CovEvent = CovEventEff), hessian = TRUE,
                fit.method = "BPM")

# Step 2: Compare informative hypotheses
# ---------------------------------------
# Load package
library(bain)						# Bayes factors for informative hypotheses

# Obtain the parameter estimates
est <- fit$coef
names(est) <- c("division", "hierarchy", "building")

# Obtain the covariance matrix
covMat <- fit$cov

# Obtain the sample size (number of events)
sampSize <- fit$m

# Define the hypotheses
hyp <- "division = hierarchy = building; division > hierarchy = building;
division > hierarchy > building; division > building > hierarchy"

# Run Bain
set.seed(8389)
results <- bain(est, hyp, n = sampSize, Sigma = covMat,
                group_parameters = 0, joint_parameters = 3)

# Print results
print(results)



fit <- rem.dyad(edgelist = relevents, n = nrow(actors), effects = "CovEvent",
                ordinal = FALSE, covar = list(CovEvent = CovEventEff), hessian = TRUE,
                fit.method = "BPM")
names(fit$coef) <- c("division", "hierarchy", "building")
BF(fit)
BF(fit,hyp)




