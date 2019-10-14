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

# Aapplication 1

# load R package BFpack and bain
library(BFpack)
library(bain)
# Frequentist one sample t test at the null point mu=5 on the therapeutic data
ttest1 <- t_test(therapeutic,mu=5)
print(ttest1)
# confirmatory Bayesian one sample t test
BF1 <- BF(ttest1,"mu=5")
summary(BF1)


# Application 2

# Fit an anova model on the tvprices data
aov1 <- aov(price ~ anchor*motivation,data=tvprices)
# Perform exploratory Bayes factor tests for ANOVA design
BF(aov1)


# Application 3
# Perform classical Bartlett test for homogeneity of group varianceson the accuracy data
bartlett <- bartlett_test(x = attention$accuracy, g = attention$group)
# Specify informative hypotheses on the group variances based on substantive expectations
hypothesis <- c("Controls=TS<ADHD;
    Controls<TS=ADHD;
    Controls=TS=ADHD")
set.seed(358)
# Perform the Bayesian hypothesis tests
BF_var <- BF(bartlett, hypothesis)
summary(BF_var)

# Application 4
# Fit the multivariate regression model on the fmri data
fmri.lm <- lm(cbind(Superficial,Middle,Deep) ~ Face + Vehicle, data=fmri)
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

# Application 4b
constraints.fmri2 <- "Face_on_Deep = Face_on_Superficial = Face_on_Middle < 0;
     Face_on_Deep < Face_on_Superficial = Face_on_Middle < 0"
fmri.lm2 <- lm(cbind(Superficial,Middle,Deep) ~ Face + Vehicle, data=fmri)
BF.fmri2 <- BF(fmri.lm2, hypothesis=constraints.fmri2)


#Application 5
# Fit logistic regression model
fit <- glm(sent ~ ztrust + zfWHR + zAfro + glasses + attract + maturity +
             tattoos, family = binomial(), data = wilson)
# Perform the Bayesian hypothesis tests
set.seed(123)
BF_glm <- BF(fit, hypothesis="ztrust > (zfWHR, zAfro) > 0; ztrust > (zfWHR, zAfro) = 0")
summary(BF_glm)


# Application 6: Correlation analysis
#Fit multivariate normal model
lm6 <- lm(cbind(Im,Del,Wmn,Cat,Fas,Rat) ~ -1 + Group, data=memory)
set.seed(123)
#Perform Bayes factor test of the order hypothesis against its compleplement
BF6_cor <- BF(lm6,parameter="correlation", hypothesis=
     "Del_with_Im_in_GroupHC > Del_with_Im_in_GroupSZ &
     Del_with_Wmn_in_GroupHC > Del_with_Wmn_in_GroupSZ &
     Del_with_Cat_in_GroupHC > Del_with_Cat_in_GroupSZ &
     Del_with_Fas_in_GroupHC > Del_with_Fas_in_GroupSZ &
     Del_with_Rat_in_GroupHC > Del_with_Rat_in_GroupSZ &
     Im_with_Wmn_in_GroupHC > Im_with_Wmn_in_GroupSZ &
     Im_with_Cat_in_GroupHC > Im_with_Cat_in_GroupSZ &
     Im_with_Fas_in_GroupHC > Im_with_Fas_in_GroupSZ &
     Im_with_Rat_in_GroupHC > Im_with_Rat_in_GroupSZ &
     Wmn_with_Cat_in_GroupHC > Wmn_with_Cat_in_GroupSZ &
     Wmn_with_Fas_in_GroupHC > Wmn_with_Fas_in_GroupSZ &
     Wmn_with_Rat_in_GroupHC > Wmn_with_Rat_in_GroupSZ &
     Cat_with_Fas_in_GroupHC > Cat_with_Fas_in_GroupSZ &
     Cat_with_Rat_in_GroupHC > Cat_with_Rat_in_GroupSZ &
     Fas_with_Rat_in_GroupHC > Fas_with_Rat_in_GroupSZ")
summary(BF6_cor)


# Application 7. Testing intraclas correlations
library(lme4)
#Only consider the timssICC data for the measurements of the four counteries in 2011
timssICC_subset <- timssICC[(timssICC$groupNL11==1)+(timssICC$groupHR11==1)+
                              (timssICC$groupDE11==1)+(timssICC$groupDK11==1)>0,]
#Fit a random intercept model with country specific random effect variances across schools
outlme1 <- lmer(math ~ -1 + gender + weight + lln +
                  groupNL11 + (0+groupNL11 | schoolID) +
                  groupHR11 + (0+groupHR11 | schoolID) +
                  groupDE11 + (0+groupDE11 | schoolID) +
                  groupDK11 + (0+groupDK11 | schoolID),
                data=timssICC_subset)
#Perform Bayes factor test on intraclass correlations across countries in 2011
set.seed(123)
BFout <- BF(outlme1,hypothesis=
     "groupNL11<groupHR11<groupDE11<groupDK11;
     groupNL11=groupHR11=groupDE11=groupDK11")
summary(BFout)
