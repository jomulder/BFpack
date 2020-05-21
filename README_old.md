# BFpack

R-functions for Bayesian exploratory (e.g., equal vs negative vs postive) and confirmatory (equality and/or order constraints) hypothesis testing under commonly used statistical models, including Bayesian t testing, Bayesian (M)AN(C)OVA, multivariate/univariate linear regression, Bayesian correlation analysis, Bayesian multilevel analysis, Bayesian generalized linear models (e.g., logistic regression). The main function `BF` needs a fitted model (e.g., an object of class `lm`) as input as well as an argument called `hypothesis` which is a string that specifies a set of equality/order constraints on the parameters. By applying the function `get_estimates`on a fitted model, the names of the parameters are returned on which constrained hypotheses can be formulated.

Developers and collaborators: Joris Mulder, Caspar van Lissa, Donald R. Williams, Xin Gu, Anton Olsson-Collentine, Florian Böing-Messing, Andrew Tomarken, Marlyne Bosman-Meijerink, Eric-Jan Wagenmakers, Yves Rosseel, Jean-Paul Fox, Janosch Menke, and Herbert Hoijtink.

Licensed under the GNU General Public License version 2 (June, 1991)


Installation
------------

You can install BFpack from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("jomulder/BFpack")


Example analyses
------------

# Application 1: Bayesian t test

# load R package BFpack and bain
install.packages("BFpack")
library(BFpack)
install.packages("bain")
library(bain)
# Frequentist one sample t test at the null point mu=5 on the therapeutic data
ttest1 <- t_test(therapeutic, alternative = "greater", mu = 5)
print(ttest1)
# confirmatory Bayesian one sample t test
BF1 <- BF(ttest1, hypothesis = "mu = 5; mu > 5", complement = FALSE)
print(BF1)
summary(BF1)


# Application 2: Bayesian analysis of variance

# Fit an anova model on the tvprices data
aov1 <- aov(price ~ anchor * motivation, data = tvprices)
# Perform exploratory Bayes factor tests for ANOVA design
BF(aov1)


# Application 3: Testing homogeneity of variances

# Perform classical Bartlett test for homogeneity of group varianceson the accuracy data
bartlett <- bartlett_test(x = attention$accuracy, g = attention$group)
# Specify informative hypotheses on the group variances based on substantive expectations
hypothesis <- c("Controls = TS < ADHD; Controls < TS = ADHD;
    Controls = TS = ADHD")
set.seed(358)
# Perform the Bayesian hypothesis tests
BF_var <- BF(bartlett, hypothesis)
summary(BF_var)


#Application 4: Bayesian logistic regression

# Fit logistic regression model
fit <- glm(sent ~ ztrust + zfWHR + zAfro + glasses + attract + maturity +
               tattoos, family = binomial(), data = wilson)
# Perform the Bayesian hypothesis tests
set.seed(123)
BF_glm <- BF(fit, hypothesis = "ztrust > (zfWHR, zAfro) > 0;
             ztrust > zfWHR = zAfro = 0")
summary(BF_glm)

ct <- lmtest::coeftest(fit)
BF(ct[,1], Sigma = diag(ct[,2]^2), n = nrow(wilson))


# Application 5: Multivariate multiple normal linear regression 

# Fit the multivariate regression model on the fmri data
fmri.lm <- lm(cbind(Superficial, Middle, Deep) ~ Face + Vehicle,
              data = fmri)
# Specify informative hypotheses on the coefficients across predictor and dependent variables
# based on substantive expectations
constraints.fmri <- "Face_on_Deep = Face_on_Superficial = Face_on_Middle < 0 <
     Vehicle_on_Deep = Vehicle_on_Superficial = Vehicle_on_Middle;
     Face_on_Deep < Face_on_Superficial = Face_on_Middle < 0 < Vehicle_on_Deep =
     Vehicle_on_Superficial = Vehicle_on_Middle"
# Perform the Bayesian hypothesis tests
set.seed(123)
BF_fmri <- BF(fmri.lm, hypothesis = constraints.fmri)
summary(BF_fmri)

# Application 4b (with random missing observations)
# Specify constrained hypotheses
constraints.fmri2 <- "Face_on_Deep = Face_on_Superficial = Face_on_Middle < 0;
     Face_on_Deep < Face_on_Superficial = Face_on_Middle < 0"
# Fit model on complete data
fmri.lm2 <- lm(cbind(Superficial,Middle,Deep) ~ Face + Vehicle, data=fmri)
# Compute Bayes factors and posterior probabilities based on complete data
BF.fmri2 <- BF(fmri.lm2, hypothesis=constraints.fmri2)

# Generate random missings
fmri_missing <- fmri
set.seed(1234)
for(i in 1:10){
    fmri_missing[sample(1:nrow(fmri), 1), sample(1:ncol(fmri), 1)] <- NA
}

# List-wise delete rows that contain at least one missing observation
fmri_listdel <- fmri_missing[!is.na(apply(fmri_missing, 1, sum)),]
fmri.lm2_listdel <- lm(cbind(Superficial, Middle, Deep) ~ Face + Vehicle,
                       data = fmri_listdel)
# Execute Bayes factor test based on list-wise deleted data
BF.fmri2_listdel <- BF(fmri.lm2_listdel, hypothesis = constraints.fmri2)
#print output
print(BF.fmri2_listdel)

# Create 500 complete impluted data sets using the mice package
M <- 500
install.packages("mice")
library(mice)
mice_fmri <- mice :: mice(data = fmri_missing, m = M, meth = c("norm",
    "norm", "norm", "norm", "norm"), diagnostics = F, printFlag = F)
# Execute tests based on complete imputed data.
# Extract the posterior and prior quantities of the extended Savage-Dickey density ratio
# from the Specification matrix for the imputed data.
relmeas_all <- matrix(unlist(lapply(1:M, function(m){
    fmri.lm_m <- lm(cbind(Superficial, Middle, Deep) ~ Face + Vehicle,
        data = mice::complete(mice_fmri, m))
    BF.fmri2_m <- BF(fmri.lm_m, hypothesis = constraints.fmri2)
    # Return four quantities of SD-ratio per constrained hypothesis
    c(BF.fmri2_m$BFtable_confirmatory[, 1:4])
    })),ncol = M)
# Compute Monte Carlo estimate of the four quantities for the constrained hypotheses
relmeas <- matrix(apply(relmeas_all, 1, mean),nrow = 3)
# Give appropriate names
row.names(relmeas) <- c("H1", "H2", "H3")
colnames(relmeas) <- c("complex=", "complex>", "fit=", "fit>")
# Compute Bayes factors using Equation 4
BF_tu_confirmatory <- relmeas[,3] * relmeas[,4] / (relmeas[,1] *
    relmeas[,2])
# Compute Posterior probabilities using Equation 3.
PHP <- BF_tu_confirmatory / sum(BF_tu_confirmatory)


# Application 6: Bayesian correlation analysis

# Create data matrices for the two separate groups.
memoryHC <- subset(memory,Group=="HC")[,-7]
memorySZ <- subset(memory,Group=="SZ")[,-7]
set.seed(123)
# Execute unconstrained Bayesian correlation analysis using cor_test
cor1 <- cor_test(memoryHC,memorySZ)
# Perform Bayes factor test of the order hypothesis against its complement
BF6_cor <- BF(cor1, hypothesis =
    "Del_with_Im_in_g1 > Del_with_Im_in_g2 &
    Del_with_Wmn_in_g1 > Del_with_Wmn_in_g2 &
    Del_with_Cat_in_g1 > Del_with_Cat_in_g2 &
    Del_with_Fas_in_g1 > Del_with_Fas_in_g2 &
    Del_with_Rat_in_g1 > Del_with_Rat_in_g2 &
    Im_with_Wmn_in_g1 > Im_with_Wmn_in_g2 &
    Im_with_Cat_in_g1 > Im_with_Cat_in_g2 &
    Im_with_Fas_in_g1 > Im_with_Fas_in_g2 &
    Im_with_Rat_in_g1 > Im_with_Rat_in_g2 &
    Wmn_with_Cat_in_g1 > Wmn_with_Cat_in_g2 &
    Wmn_with_Fas_in_g1 > Wmn_with_Fas_in_g2 &
    Wmn_with_Rat_in_g1 > Wmn_with_Rat_in_g2 &
    Cat_with_Fas_in_g1 > Cat_with_Fas_in_g2 &
    Cat_with_Rat_in_g1 > Cat_with_Rat_in_g2 &
    Fas_with_Rat_in_g1 > Fas_with_Rat_in_g2")
print(BF6_cor)


# Application 7. Bayesian testing of intraclas correlations
#Only consider the timssICC data for the measurements of the four counteries in 2011
timssICC_subset <- subset(timssICC, groupNL11 == 1 | groupHR11 == 1 |
    groupDE11 == 1 | groupDK11 == 1)
#Fit a random intercept model with country specific random effect variances across schools
outlme1 <- lme4::lmer(math ~ -1 + gender + weight + lln +
                  groupNL11 + (0+groupNL11 | schoolID) +
                  groupHR11 + (0+groupHR11 | schoolID) +
                  groupDE11 + (0+groupDE11 | schoolID) +
                  groupDK11 + (0+groupDK11 | schoolID),
                data=timssICC_subset)
#Perform Bayes factor test on intraclass correlations across countries in 2011
set.seed(123)
BFicc <- BF(outlme1,hypothesis=
     "groupNL11<groupHR11<groupDE11<groupDK11;
     groupNL11=groupHR11=groupDE11=groupDK11")
summary(BFicc)
round(BFicc$estimates,3)

