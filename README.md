# BFpack

R-functions for Bayesian model selection and hypothesis testing using default Bayes factors. Currently, the package only contains a functions for computing the BIC for order constrianed models (Mulder & Raftery, in prep.). The functions needs a fitted model (e.g., glm, survival) as input as well as a string that specifies a set of order constraints on the regression coefficients.


Basic example
-------------

``` r
library(BFpack)

#An model selection problem between a model assuming a positive effect, a negative effect, or no effect.

salary <- read.table("http://data.princeton.edu/wws509/datasets/salary.dat", header=TRUE)
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


```

Installation
------------

You can install BFpack from github with:

``` r
# install.packages("devtools")
devtools::install_github("jomulder/BFpack")
