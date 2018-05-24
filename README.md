# BFpack

R-functions for Bayesian model selection and hypothesis testing using default Bayes factors. Currently, the package only contains a functions for computing the BIC for order constrianed models (Mulder & Raftery, in prep.). The functions needs a fitted model (e.g., glm, survival) as input as well as a string that specifies a set of order constraints on the regression coefficients.


Basic example
-------------

``` r
library(BFpack)


salary <- read.table("http://data.princeton.edu/wws509/datasets/salary.dat", header=TRUE)

# Testing a model assuming a positive gender effect versus a model assuming a negative
# gender effect versus a model assuming no gender effect.

#fit a model with all effects present
salfit <- glm( sl ~ sx + rk + yr + dg + yd, family = gaussian, data=salary)
#a model which assumes a positive effect for being male
bic_oc(salfit,"sxmale>0")
#a model which assumes a negative effect for being male
bic_oc(salfit,"sxmale<0")

#fit a model excluding the gender effect
salfit0 <- glm( sl ~ rk + yr + dg + yd, family = gaussian, data=salary)
#a model which assumes no gender effect
bic_oc(salfit0)

#getting posterior probabilities of the models assuming they are equally likely a priori
bicvec <- c(bic_oc(salfit,"sxmale>0"),bic_oc(salfit,"sxmale<0"),bic_oc(salfit0))
postprob(bicvec)


#testing an ordered effect that an associate professor has a higher salary than an assistent professor, and the a full professor has a higher salary than an associate professor in this population
bic_oc(salfit,"rkfull>rkassociate>0")
bic_oc(salfit,"rkfull>rkassociate>0",complement=TRUE)
```

Installation
------------

You can install BFpack from github with:

``` r
# install.packages("devtools")
devtools::install_github("jomulder/BFpack")
