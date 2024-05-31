#' BFpack: Flexible Bayes factor testing of scientific expectations
#'
#'
#' @description The \code{R} package \strong{BFpack} provides tools for exploratory and
#' confirmatory Bayesian hypothesis testing using Bayes factors and posterior probabilities
#' under common statistical models. The main function `BF` needs a fitted model `x` as input
#' argument. Depending on the class of the fitted model, a standard hypothesis test is
#' executed by default. For example, if `x` is a
#' fitted regression model of class `lm` then posterior probabilities are computed of whether
#' each separate coefficient is zero, negative, or positive (assuming equal prior probabilities).
#' If one has specific hypotheses with
#' equality and/or order constraints on the parameters under the fitted model `x` then these
#' can be formulated using the `hypothesis` argument (a character string), possibly together
#' prior probabilities for the hypotheses via the `prior` argument (default all hypotheses are
#' equally likely a priori), and the `complement` argument which is a logical stating whether
#' the complement hypotheses should be included in the case (`TRUE` by default).
#'
#' Use compilation for Fortran functions
#'
#' @references
#' Mulder, J., D.R. Williams, Gu, X., A. Tomarken,
#' F. Böing-Messing, J.A.O.C. Olsson-Collentine, Marlyne Meyerink, J. Menke,
#' J.-P. Fox, Y. Rosseel, E.J. Wagenmakers, H. Hoijtink., and van Lissa, C.
#' (submitted). BFpack: Flexible Bayes Factor Testing of Scientific Theories
#' in R. \url{https://arxiv.org/abs/1911.07728}
#'
#' Mulder, J., van Lissa, C., Gu, X., Olsson-Collentine, A., Boeing-Messing, F., Williams,
#' D. R., Fox, J.-P., Menke, J., et al. (2019). BFpack: Flexible Bayes Factor Testing of
#' Scientific Expectations. (Version 0.2.1) \url{https://CRAN.R-project.org/package=BFpack}
#'
#'
#' @examples
#' \dontrun{
#' # EXAMPLE 1. One-sample t test
#' ttest1 <- t_test(therapeutic, mu = 5)
#' print(ttest1)
#' # confirmatory Bayesian one sample t test
#' BF1 <- BF(ttest1, hypothesis = "mu = 5")
#' summary(BF1)
#' # exploratory Bayesian one sample t test
#' BF(ttest1)
#'
#' # EXAMPLE 2. ANOVA
#' aov1 <- aov(price ~ anchor * motivation,data = tvprices)
#' BF1 <- BF(aov1, hypothesis = "anchorrounded = motivationlow;
#'                               anchorrounded < motivationlow")
#' summary(BF1)
#'
#' # EXAMPLE 3. Logistic regression
#' fit <- glm(sent ~ ztrust + zfWHR + zAfro + glasses + attract + maturity +
#'    tattoos, family = binomial(), data = wilson)
#' BF1 <- BF(fit, hypothesis = "ztrust > zfWHR > 0;
#'                              ztrust > 0 & zfWHR = 0")
#' summary(BF1)
#'
#' # EXAMPLE 4. Correlation analysis
#' set.seed(123)
#' cor1 <- cor_test(memory[1:20,1:3])
#' BF1 <- BF(cor1)
#' summary(BF1)
#' BF2 <- BF(cor1, hypothesis = "Wmn_with_Im > Wmn_with_Del > 0;
#'                               Wmn_with_Im = Wmn_with_Del = 0")
#' summary(BF2)
#' }
#'
#' @docType package
#'
#'
#' @name BFpack-package
#'
"_PACKAGE"
