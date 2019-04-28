# BFpack

R-functions for Bayesian exploratory and confirmatory hypothesis testing for the most commonly used statistical models, including (but not limited to) univariate/multivariate t testing, (M)AN(C)OVA, multivariate/univariate regression, structural equation modeling, (mixed) generalized linear models. The functions need fitted models (e.g., lm) as input as well as a string that specifies a set of order constraints on the regression coefficients.

Developers and collaborators: Joris Mulder, Caspar van Lissa, Xin Gu, Anton Olsson-Colentine, Florian Böing-Messing, Andrew Tomarken, Eric-Jan Wagenmakers, Herbert Hoitink.

Licensed under the GNU General Public License version 2 (June, 1991)


Installation
------------

You can install BFpack from github with:

``` r
# install.packages("devtools")
devtools::install_github("jomulder/BFpack")
