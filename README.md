# BFpack

R-functions for Bayesian exploratory (equal vs negative vs postive) and confirmatory (equality and/or order constraints) hypothesis testing for commonly used statistical models, including (but not limited to) univariate/multivariate t testing, (M)AN(C)OVA, multivariate/univariate regression, random intercept models. The functions need fitted models (e.g., lm) as input as well as a string that specifies a set of order constraints on the regression coefficients.

Developers and collaborators: Joris Mulder, Caspar van Lissa, Xin Gu, Anton Olsson-Collentine, Florian Böing-Messing, Donald R. Williams, Andrew Tomarken, Marlyne Bosman-Meijerink, Eric-Jan Wagenmakers, Yves Rosseel, Jean-Paul Fox, Janosch Menke, and Herbert Hoijtink.

Licensed under the GNU General Public License version 2 (June, 1991)


Installation
------------

You can install BFpack from github with:

``` r
# install.packages("devtools")
devtools::install_github("jomulder/BFpack")


Example analyses
------------

# load R package BFpack and bain
library(BFpack)
library(bain)

# EXAMPLE 1. One-sample t test
ttest1 <- t_test(therapeutic,mu=5)
print(ttest1)
# confirmatory Bayesian one sample t test
BF1 <- BF(ttest1,"mu=5")
summary(BF1)
# exploratory Bayesian one sample t test
BF(ttest1)

# EXAMPLE 2. ANOVA
aov1 <- aov(price ~ anchor*motivation,data=tvprices)
BF1 <- BF(aov1,hypothesis="anchorrounded=motivationlow;
   anchorrounded<motivationlow;
   anchorrounded>motivationlow")
   summary(BF1)

# EXAMPLE 3. Logistic regression
fit <- glm(sent ~ ztrust + zfWHR + zAfro + glasses + attract + maturity +
   tattoos, family = binomial(), data = wilson)
BF1 <- BF(fit, hypothesis = "ztrust > (zfWHR, zAfro) > 0;
    ztrust > (zfWHR, zAfro) = 0")
summary(BF1)

# EXAMPLE 4. Correlation analysis
res <- polycor::hetcor(fmri[,3:5])
BF1 <- BF(res)
summary(BF1)
BF1 <- BF(res,hypothesis="(Middle_with_Superficial,Deep_with_Superficial,Deep_with_Middle) > 0;
          Middle_with_Superficial=Deep_with_Superficial=Deep_with_Middle=0")
summary(BF1)





