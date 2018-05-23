# BFpack

R-functions for Bayesian model selection and hypothesis testing using default Bayes factors.

Currently, the package only contains a functions for computing the BIC for order constrianed models (Mulder & Raftery, in prep.). The functions needs a fitted model (e.g., glm, survival) as input as well as a string that specifies a set of order constraints on the regression coefficients.


Basic example
-------------

``` r
library(BFpack)

EXAMPLE HERE


```

Installation
------------

You can install BFpack from github with:

``` r
# install.packages("devtools")
devtools::install_github("jomulder/BFpack")
