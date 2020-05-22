
<img src="man/figures/logo_BFpack.png" width = 850 />

# 

[![CRAN
Version](http://www.r-pkg.org/badges/version/BFpack)](https://cran.r-project.org/package=BFpack)
[![Downloads](https://cranlogs.r-pkg.org/badges/BGGM)](https://cran.r-project.org/package=BFpack)

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

Below several example analyses are provided using **BFpack**.

### Bayesian t testing

First a classical one sample t test is executed for the test value
\(\mu = 5\) on the therapeutic data

``` r
ttest1 <- bain::t_test(therapeutic, alternative = "greater", mu = 5)
```

The `t_test` function is part of the ***bain*** package. The function is
equivalent to the standard `t.test` function with the addition that the
returned object contains additional output than the standard `t.test`
function.

To perform a Bayesian t test plug the fitted object into the `BF`
function.

``` r
library(BFpack)
BF1 <- BF(ttest1)
```

This executes an exhaustive test around the null value: `H1: mu = 5`
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

### Correlation analysis

By default `BF` performs exhaustice tests of whether the separate
correlations are zero, negative, or positive. The name of the
correlations is constructed using the names of the variables separated
by `_with_`.

``` r
set.seed(123)
cor1 <- cor_test(memory[,1:3])
BF1 <- BF(cor1)
print(BF1)
```

Constraints can also be tested between correlations, e.g., whether all
correlations are equal and positive versus an unconstrained
complement.

``` r
BF2 <- BF(cor1, hypothesis = "Del_with_Im = Wmn_with_Im = Wmn_with_Del > 0")
print(BF2)
```

### Univariate/Multivariate multiple regression

For a univariate regression model, by default an exhaustive test is
executed of whether an effect is zero, negative, or postive.

``` r
lm1 <- lm(Superficial ~ Face + Vehicle, data = fmri)
BF1 <- BF(lm1)
print(BF1)
```

Hypotheses can be tested with equality and/or order constraints on the
effects of interest. If prefered the complement hypothesis can be
omitted using the `complement`
argument

``` r
BF2 <- BF(lm1, hypothesis = "Vehicle > 0 & Face < 0; Vehicle = Face = 0",
          complement = FALSE)
print(BF2)
```

In a multivariate regression model hypotheses can be tested on the
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

## Citing BFpack

You can cite the package and the paper using the following reference

> Mulder, J., van Lissa, C., Gu, X., Olsson-Collentine, A.,
> Boeing-Messing, F., Williams, D. R., Fox, J.-P., Menke, J., et al.
> (2019). BFpack: Flexible Bayes Factor Testing of Scientific
> Expectations. (Version 0.2.1) \[R package\].
> <https://CRAN.R-project.org/package=BFpack>

> Mulder, J., Gu, X., Olsson-Collentine, A., Tomarken, A.,
> Böing-Messing, F., Hoijtink, H., . . . van Lissa, C. (2019). BFpack:
> Flexible Bayes factor testing of scientific theories in R. Retrieved
> from <https://arxiv.org/abs/1911.07728>

## Contributing and Contact Information

If you have ideas, please get involved. You can contribute by opening an
issue on GitHub, or sending a pull request with proposed features.

  - File a GitHub issue [here](https://github.com/jomulder/BFpack)
  - Make a pull request [here](https://github.com/jomulder/BFpack/pulls)

By participating in this project, you agree to abide by the [Contributor
Code of Conduct v2.0](code_of_conduct.md).
