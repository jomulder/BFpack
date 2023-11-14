
<img src="man/figures/logo_BFpack.png" width = 850 />

# 

[![CRAN
Version](http://www.r-pkg.org/badges/version/BFpack)](https://cran.r-project.org/package=BFpack)
[![Downloads](https://cranlogs.r-pkg.org/badges/BFpack)](https://cran.r-project.org/package=BFpack)
[![R-CMD-check](https://github.com/jomulder/BFpack/workflows/R-CMD-check/badge.svg)](https://github.com/jomulder/BFpack/actions)
<!-- Insert codecov badge here -->

[![Contributor
Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](https://www.contributor-covenant.org/version/2/0/code_of_conduct.html)

The `R` package **BFpack** contains a set of functions for exploratory
hypothesis testing (e.g., equal vs negative vs postive) and confirmatory
hypothesis testing (with equality and/or order constraints) using Bayes
factors and posterior probabilities under commonly used statistical
models, including (but not limited to) Bayesian t testing, (M)AN(C)OVA,
multivariate/univariate linear regression, correlation analysis,
multilevel analysis, or generalized linear models (e.g., logistic
regression). The main function `BF` needs a fitted model (e.g., an
object of class `lm` for a linear regression model) and (optionally) the
argument `hypothesis`, a string which specifies a set of equality/order
constraints on the parameters. By applying the function
`get_estimates`on a fitted model, the names of the parameters are
returned on which constrained hypotheses can be formulated. Bayes
factors and posterior probabilities are computed for the hypotheses of
interest.

## Installation

Install the latest release version of `BFpack` from CRAN:

``` r
install.packages("BFpack")
```

The current developmental version can be installed with

``` r
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   
remotes::install_github("jomulder/BFpack")
```

## Example analyses

Below several example analyses are provided using **BFpack**. As input the main function `BF()` requires
a fitted model which which the necessary elements are extracted to compute Bayes factors and posterior
probabilities for the hypotheses.

### Bayesian t testing

#### Univariate t testing

First a classical one sample t test needs executed on the test value
\(mu = 5\) on the `therapeutic` data (part of `BFpack`). Here a right one-tailed classical test is executed:

``` r
ttest1 <- bain::t_test(therapeutic, alternative = "greater", mu = 5)
```

The `t_test` function is part of the **bain** package. The function is
equivalent to the standard `t.test` function with the addition that the
returned object contains additional output than the standard `t.test`
function.

Two default Bayes factors are implemented in **BFpack** to execute a t test: the fractional Bayes factor
(O'Hagan, 1995) and the prior adjusted fractional Bayes factor (Mulder, 2014). Both do not require a
prior to be manually specified as a default prior is implicitly constructed using a minimal fraction
of the data. The remaining fraction is used for hypothesis testing. The fractional Bayes factor
behaves similar as the JZS Bayes factor (Rouder et al., 2009, as implemented in the **BayesFactor** package)
for standard null hypothesis testing and the prior adjusted Bayes factor was specifically designed for
testing one-sided hypotheses. When setting the argument `BF.type=1` or `BF.type=2`, the fractional
Bayes factor and the prior adjusted fractional Bayes factor is used, respectively. The default choice
is `BF.type=2`.

To perform a Bayesian t test, the `BF` function is run on the fitted object. 

``` r
library(BFpack)
BF1 <- BF(ttest1)
```

This executes an exploratory test around the null value: `H1: mu = 5`
versus `H2: mu < 5` versus `H3: mu > 5` assuming equal prior
probabilities for `H1`, `H2`, and `H3` of 1/3. The output presents the
posterior probabilities for the three hypotheses.

The same test would be executed when the same hypotheses are explicitly
specified using the `hypothesis` argument.

``` r
hypothesis <- "mu = 5; mu < 5; mu > 5"
BF(ttest1, hypothesis = hypothesis)
```

When testing hypotheses via the `hypothesis` argument, the output also
presents an `Evidence matrix` containing the Bayes factors between the
hypotheses.

The argument `prior.hyp` can be used to specify different prior probabilities
for the hypotheses. For example, when the left one-tailed hypothesis is not possible
based on prior considerations (e.g., see [preprint](https://arxiv.org/abs/1911.07728)) while the precise (null) hypothesis and the right
one-tailed hypothesis are equally likely, the argument `prior.hyp` should be a vector
specifying the prior probabilities of the respective hypotheses
``` r
BF(ttest1, hypothesis = "mu = 5; mu < 5; mu > 5", prior.hyp = c(.5,0,.5))
```

#### Multivariate t testing

Bayesian multivariate t tests can be executed by first fitting a multivariate (regression) model using the `lm` function, and subsequently, the means of the dependent variables (or other coefficients) in the model can be tested using the `BF()` function. Users have to be aware that the (adjusted) means are modeled using intercepts which are named `(Intercept)` by default by `lm` while the `hypothesis` argument in `BF()` does not allow effect names that include brackets (i.e., `(` or `)`). To circumvent this, one can create a vector of 1s, with name (say) `ones`, which can replace the intercept but results in an equivalent model.

For example, let us consider a multivariate normal model for the dependent variables `Superficial`, `Middle`, and `Deep` in the `fmri` data set:

``` r
fmri1 <- cbind(fmri,ones=1)
mlm1 <- lm(cbind(Superficial,Middle,Deep) ~ -1 + ones, data = fmri1)
```

Next, we can (for instance) test whether all means are equal to 0 (`H1`), whether all means are positive (`H2`), or neither (`complement`):

``` r
BFmlm1 <- BF(mlm1,
             hypothesis="ones_on_Superficial=ones_on_Middle=ones_on_Deep=0;
                        (ones_on_Superficial,ones_on_Middle,ones_on_Deep)>0",
             complement = TRUE)
```

### Analysis of variance

First an analysis of variance (ANOVA) model is fitted using the `aov`
fuction in `R`.

``` r
aov1 <- aov(price ~ anchor * motivation, data = tvprices)
```

Next a Bayesian test can be performed on the fitted object.

``` r
BF(aov1)
```

By default posterior probabilities are computed of whether main effects
and interaction effects are present. Alternative constrained hypotheses
can be tested on the model parameters `get_estimates(aov1)`.

Similar as for the Bayesian t test, two default Bayes factors are implemented in **BFpack** for (multivariate)
analysis of (co)variance: the fractional Bayes factor (O'Hagan, 1995) and the prior adjusted fractional Bayes factor (Mulder, 2014, Mulder and Gu, 2022). The default choice is `BF.type=2`, which uses the prior adjusted fractional
Bayes factor.


### Univariate/Multivariate multiple regression

As an example we consider the `fmri` data (McGuigin et al, 2020) as discussed in Mulder et al. (2021). First, a classical linear regression is fitted with dependent variable `Deep` and predictor variables `Face` and `Vehicle`:

``` r
lm1 <- lm(Deep ~ Face + Vehicle, data = fmri)
```

When not using the `hypothesis` argument, Bayes factors and posterior probabilities are given of whether each predictor
has a zero, negative, or positive effect (assuming equal prior probabilities) against the full model:

``` r
BF(lm1)
```

As motivated in Mulder et al. (2021), it was expected that `Face` had a negative effect on `Deep` and `Vehicle` had a positive effect on `Deep`. This (combined) one-sided hypothesis can be tested against its complement according to

``` r
BF(lm1, hypothesis = "Face < 0 < Vehicle")
```

The hypothesis of interest receives clear support from the data.

In a multivariate multiple regression model, hypotheses can be tested on the
effects on the same dependent variable, and on effects across different
dependent variables. The name of an effect is constructed as the name of
the predictor variable and the dependent variable separated by `_on_`.
Testing hypotheses with both constraints within a dependent variable and
across dependent variables makes use of a Monte Carlo estimate which may
take a few seconds.

``` r
lm2 <- lm(cbind(Superficial, Middle, Deep) ~ Face + Vehicle,
              data = fmri)
constraint2 <- "Face_on_Deep = Face_on_Superficial = Face_on_Middle < 0 <
     Vehicle_on_Deep = Vehicle_on_Superficial = Vehicle_on_Middle;
     Face_on_Deep < Face_on_Superficial = Face_on_Middle < 0 < Vehicle_on_Deep =
     Vehicle_on_Superficial = Vehicle_on_Middle"
set.seed(123)
BF3 <- BF(lm2, hypothesis = constraint2)
summary(BF3)
```

Finally note that for (multivariate) multiple regression again note that two default Bayes factors are implemented in **BFpack**: the fractional Bayes factor (O'Hagan, 1995) and the prior adjusted fractional Bayes factor (Mulder, 2014; Mulder and Gu, 2022) which can be chosen using the argument `BF.type=1` and `BF.type=2`, respectively. The default
choice is the prior adjusted fractional Bayes factor. This criterion was specifically designed for testing one-sided
and order constrained hypotheses.


### Logistic regression

An example hypothesis test is consdered under a logistic regression
model. First a logistic regression model is fitted using the `glm`
function

``` r
fit_glm <- glm(sent ~ ztrust + zfWHR + zAfro + glasses + attract + maturity +
               tattoos, family = binomial(), data = wilson)
```

The names of the regression coefficients on which constrained hypotheses
can be formualted can be extracted using the `get_estimates` function.

``` r
get_estimates(fit_glm)
```

Two different hypotheses are formulated with competing equality and/or
order constraints on the parameters of interest. These hypotheses are
motivated in Mulder et al. (2019)

``` r
BF_glm <- BF(fit_glm, hypothesis = "ztrust > (zfWHR, zAfro) > 0;
             ztrust > zfWHR = zAfro = 0")
summary(BF_glm)
```

By calling the `summary` function on the output object of class `BF`,
the results of the exploratory tests are presented of whether each
separate parameter is zero, negative, or positive, and the results of
the confirmatory test of the hypotheses under the `hypothesis` argument
are presented. When the hypotheses do not cover the complete parameter
space, by default the complement hypothesis is added which covers the
remaining parameter space that is not covered by the constraints under
the hypotheses of interest. In the above example, the complement
hypothesis covers the parameter space where neither `"ztrust > (zfWHR,
zAfro) > 0"` holds, nor where `"ztrust > zfWHR = zAfro = 0"` holds.

The Bayes factors and posterior posterior probabilities are based on the
approximated adjusted default Bayes factor (Gu et al., 2018).


### Correlation analysis

Bayes factors and posterior posterior probabilities among constrained hypotheses on measures of association
are computed using uniform prior for the correlations (Mulder and Gelissen, 2023).
By default `BF` performs exhaustice tests of whether the separate
correlations are zero, negative, or positive. The name of the
correlations is constructed using the names of the variables separated
by `_with_`. To compute Bayes factors and posterior probabilities, first
the unconstrained model needs to be fit using the `cor_test()` function. The
resulting object can be added to the `BF()` function:

``` r
set.seed(123)
cor1 <- cor_test(memory[,1:3])
BF1 <- BF(cor1)
print(BF1)
```

Constraints can also be tested between correlations, e.g., whether all
correlations are equal and positive versus an unconstrained complement. The function
`get_estimates()` gives the names of the correlations on which constrained hypotheses
can be formulated:

``` r
get_estimates(cor1)
BF2 <- BF(cor1, hypothesis = "Del_with_Im = Wmn_with_Im = Wmn_with_Del > 0")
print(BF2)
```

Depending on the class of the variables (`numeric`, `ordered`, `factor` with 2 levels), biserial,
polyserial, polychoric, tetrachoric, or product-moment correlations are tested. As an illustration
of these other types of measures of association, we change the measurement levels of a subset of
the `mtcars` data:

``` r
mtcars_test <- mtcars[,c(1,2,9,10)]
mtcars_test[,2] <- as.ordered(mtcars_test[,2])
mtcars_test[,3] <- as.factor(mtcars_test[,3])
mtcars_test[,4] <- as.integer(mtcars_test[,4])
```

To compute the Bayes factors and posterior probabilities, again we first the full unconstrained model
using the `cor_test()` function. The resulting object is placed in the `BF()` function to obtain the
Bayes factors and posterior proabilities:

``` r
cor2 <- cor_test(mtcars_test)
BF2 <- BF(cor2, hypothesis = "0 < am_with_mpg = gear_with_mpg")
print(BF2)
```


### Running `BF` on a named vector

The input for the `BF` function can also be a named vector containing
the estimates of the parameters of interest. In this case the error
covariance matrix of the estimates is also needed via the `Sigma`
argument, as well as the sample size that was used for obtaining the
estimates via the `n` argument. Bayes factors are then computed using
Gaussian approximations of the likelihood (and posterior), similar as in
classical Wald test.

We illustrate this for a Poisson regression model

``` r
poisson1 <- glm(formula = breaks ~ wool + tension, data = datasets::warpbreaks,
             family = poisson)
```

The estimates, the error covariance matrix, and the sample size are
extracted from the fitted model

``` r
estimates <- poisson1$coefficients
covmatrix <- vcov(poisson1)
samplesize <- nobs(poisson1)
```

Constrained hypotheses on the parameters `names(estimates)` can then be
tested as follows

``` r
BF1 <- BF(estimates, Sigma = covmatrix, n = samplesize, hypothesis = 
  "woolB > tensionM > tensionH; woolB = tensionM = tensionH")
```

Note that the same hypothesis test would be executed when calling

``` r
BF2 <- BF(poisson1, hypothesis = "woolB > tensionM > tensionH;
          woolB = tensionM = tensionH")
```

because the same Bayes factor is used when running `BF` on an object of
class `glm` (see `Method: Bayes factor using Gaussian approximations`
when calling `print(BF11)` and `print(BF2)`).

The Bayes factors and posterior posterior probabilities on named vectors
are based on the adjusted default Bayes factor using Gaussian approximations
(Gu et al., 2018).


## Citing **BFpack**

You can cite the package and the paper using the following reference

**Software paper**
> Mulder, J., Williams, D. R., Gu, X., Olsson-Collentine, A., Tomarken,
> A., Böing-Messing, F., Hoijtink, H., . . . van Lissa  (2021).
> BFpack: Flexible Bayes factor testing of scientific theories in R.
> Journal of Statistical Software, 100(18).
> Retrieved from <https://arxiv.org/abs/1911.07728>

**Software package**
> Mulder, J., van Lissa, C., Gu, X., Olsson-Collentine, A.,
> Boeing-Messing, F., Williams, D. R., Fox, J.-P., Menke, J., et
> al. (2020). BFpack: Flexible Bayes Factor Testing of Scientific
> Expectations. (Version 1.3.0) \[R package\].
> <https://CRAN.R-project.org/package=BFpack>

Other references with technical details of the methodology, please see

**Bayes factors for (multivariate) t tests, (M)AN(C)OVA, (multivariate) regression**
> Mulder, J. (2014). Prior adjusted default Bayes factors for testing (in)equality
> constrained hypotheses. Computational Statistics and Data Analysis, 71, 448–463.
> <https://doi.org/10.1016/j.csda.2013.07.017>

> Mulder, J. & Gu, X. (2022) Bayesian Testing of Scientific Expectations
> under Multivariate Normal Linear Models, Multivariate Behavioral Research,
57:5, 767-783.
> <https://doi.org/10.1080/00273171.2021.1904809>

**Bayes factors for testing measures of association (e.g., correlations)** 
> Mulder, J., & Gelissen, J. P. (2023). Bayes factor testing of equality
> and order constraints on measures of association in social research.
> Journal of Applied Statistics, 50(2), 315-351.
> <https://doi.org/10.1080/02664763.2021.1992360>

**Default Bayes factors using Gaussian approximations**
> Gu. X., Mulder, J., & Hoijtink, J. (2018). Approximated adjusted
> fractional Bayes factors: A general method for testing informative
> hypotheses. British Journal of Mathematical and Statistical
> Psychology.
> <https://doi.org/10.1111/bmsp.12110>

**Bayes factors under exponential random graphs**
> Mulder, J., Friel, N., & Leifeld, P. (2023). Bayesian Testing of Scientific
> Expectations Under Exponential Random Graph Models.
> <https://doi.org/10.48550/arXiv.2304.14750>

**Bayes factors of intraclass correlations**
> Mulder, J., & Fox, J.-P. (2019). Bayes Factor Testing of Multiple
> Intraclass Correlations. Bayesian Analysis. 14(2), 521-552.
> <https://doi.org/10.1214/18-BA1115>

**Bayes factors for testing group variances**
> Böing-Messing, F., van Assen, M. A. L. M., Hofman, A. D., Hoijtink, H., 
> & Mulder, J. (2017). Bayesian evaluation of constrained hypotheses on
> variances of multiple independent groups. Psychological Methods, 22(2), 
> 262–287.
> <https://doi.org/10.1037/met0000116>

> Böing-Messing, F. & Mulder, J. (2018). Automatic Bayes factors for testing
> equality-and inequality-constrained hypotheses on variances. Psychometrika,
> 83, 586–617.
> <https://link.springer.com/article/10.1007/s11336-018-9615-z>

**Bayes factors for meta-analyses**
> Van Aert R.C.M. & Mulder, J. (2022). Bayesian hypothesis testing and 
> estimation under the marginalized random-effects meta-analysis model.
> Psychonomic Bulletin & Review, 29, 55–69.
> <https://doi.org/10.3758/s13423-021-01918-9>


## Contributing and Contact Information

If you have ideas, please get involved. You can contribute by opening an
issue on GitHub, or sending a pull request with proposed features.

  - File a GitHub issue [here](https://github.com/jomulder/BFpack)
  - Make a pull request [here](https://github.com/jomulder/BFpack/pulls)

By participating in this project, you agree to abide by the Contributor
Code of Conduct v2.0.
