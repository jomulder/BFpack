#' @title Bayes factors for Bayesian exploratory and confirmatory hypothesis
#' testing
#' @description The \code{BF} function can be used for hypothesis testing and
#'  model
#' selection using the Bayes factor. By default exploratory hypothesis tests are
#' performed of whether each model parameter equals zero, is negative, or is
#'  positive.
#' Confirmatory hypothesis tests can be executed by specifying hypotheses with
#' equality and/or order constraints on the parameters of interest.
#'
#' @param x An R object containing the outcome of a statistical analysis.
#' @param hypothesis A character string containing the informative hypotheses to
#' evaluate. The default is NULL, which will result in an exploratory analysis.
#' @param prior A vector specifying the prior probabilities of the hypotheses.
#' The default is NULL which will specify equal prior probabilities.
#' @param parameter A characted string specifying which type of parameter is of
#' interest. For a fitted model of class \code{mlm}, this can be either
#' "regression" or "correlation", where the first is the default setting.
#' @param ... Parameters passed to and from other functions.
#' @return The output is an object of class \code{BF}. The object has elements:
#' BFtu_exploratory, PHP_exploratory, BFtu_confirmatory, PHP_confirmatory,
#' BFmatrix_confirmatory, BFtable_confirmatory, BFtu_main, PHP_main,
#' BFtu_interaction, PHP_interaction, prior, hypotheses, estimates, model, call.
#' @details The function requires a fitted modeling object. Current analyses
#' that are supported: \code{\link[bain]{t_test}},
#' \code{\link[BFpack]{bartlett_test}},
#' \code{\link[stats]{aov}}, \code{\link[stats]{manova}},
#' \code{\link[stats]{lm}}, \code{mlm},
#' \code{\link[stats]{glm}}, \code{\link[polycor]{hetcor}},
#' \code{\link[lme4]{lmer}}, \code{\link[survival]{coxph}},
#' \code{\link[survival]{survreg}},
#' \code{\link[pscl]{zeroinfl}}, and \code{\link[MASS]{polr}}.
#'
#' For testing means and regression coefficients of model classes \code{t_test},
#' \code{aov}, \code{manova}, \code{lm}, \code{mlm}, adjusted fractional Bayes
#' factors
#' are computed. For testing group variances using
#' \code{\link[BFpack]{bartlett_test}},
#' adjusted fractional Bayes factors are computed.
#' For testing measures of association (e.g., correlations) under
#' model class \code{mlm} and for testing intraclass correlations under model
#' class
#' \code{lmerMod}, default Bayes factors based on uniform priors are computed.
#' For
#' all other model classes an approximate Bayes factor is computed using a
#' Gaussian
#' approximation of the posterior, similar as a classical Wald test.
#'
#' @references Mulder, J., Gu, X., A. Tomarken, F. Boing-Messing,
#' J.A.O.C. Olsson-Collentine, Marlyne Meyerink, D.R. Williams, J. Menke,
#' J.-P. Fox, Y. Rosseel, E.J. Wagenmakers, H. Hoijtink., and van Lissa, C.
#' (submitted). BFpack: Flexible Bayes Factor Testing of Scientific Theories
#' in R.
#' @examples
#' \dontshow{
#' # EXAMPLE 1. One-sample t test
#' ttest1 <- t_test(therapeutic,mu=5)
#' print(ttest1)
#' # confirmatory Bayesian one sample t test
#' BF1 <- BF(ttest1,"mu=5")
#' summary(BF1)
#' # exploratory Bayesian one sample t test
#' BF(ttest1)
#'
#' # EXAMPLE 2. ANOVA
#' aov1 <- aov(price ~ anchor*motivation,data=tvprices)
#' BF1 <- BF(aov1,hypothesis="anchorrounded=motivationlow;
#'    anchorrounded<motivationlow")
#' summary(BF1)
#'
#'
#' # EXAMPLE 3. Logistic regression
#' fit <- glm(sent ~ ztrust + zfWHR + zAfro + glasses + attract + maturity +
#'    tattoos, family = binomial(), data = wilson)
#' BF1 <- BF(fit, hypothesis = "ztrust > zfWHR > 0;
#'                              ztrust > 0 & zfWHR = 0")
#' summary(BF1)
#'
#' # EXAMPLE 4. Correlation analysis
#' res <- polycor::hetcor(fmri[,3:4])
#' BF1 <- BF(res)
#' summary(BF1)
#' BF1 <- BF(res,hypothesis="Middle_with_Superficial > 0;
#'                           Middle_with_Superficial=0")
#' summary(BF1)
#' }
#' \donttest{
#' # EXAMPLE 1. One-sample t test
#' ttest1 <- bain::t_test(therapeutic,mu=5)
#' print(ttest1)
#' # confirmatory Bayesian one sample t test
#' BF1 <- BF(ttest1,"mu=5")
#' summary(BF1)
#' # exploratory Bayesian one sample t test
#' BF(ttest1)
#'
#'
#' # EXAMPLE 2. ANOVA
#' aov1 <- aov(price ~ anchor*motivation,data=tvprices)
#' # check the names of the model parameters
#' names(aov1$coefficients)
#' BF1 <- BF(aov1,hypothesis="anchorrounded=motivationlow;
#'                            anchorrounded<motivationlow;
#'                            anchorrounded>motivationlow")
#' summary(BF1)
#'
#'
#' # EXAMPLE 3. Logistic regression
#' fit <- glm(sent ~ ztrust + zfWHR + zAfro + glasses + attract + maturity +
#'    tattoos, family = binomial(), data = wilson)
#' BF1 <- BF(fit, hypothesis = "ztrust > (zfWHR, zAfro) > 0;
#'                              ztrust > 0 & zfWHR=zAfro= 0")
#' summary(BF1)
#'
#' # EXAMPLE 4. Correlation analysis
#' res <- polycor::hetcor(fmri[,3:5])
#' BF1 <- BF(res)
#' summary(BF1)
#' BF1 <- BF(res,hypothesis="(Middle_with_Superficial,Deep_with_Superficial,
#'                           Deep_with_Middle) > 0;
#'                           Middle_with_Superficial=Deep_with_Superficial=
#'                           Deep_with_Middle=0")
#' summary(BF1)
#' }
#' @rdname BF
#' @export
#' @useDynLib BFpack, .registration = TRUE
#'
BF <- function(x, hypothesis, prior, parameter, ...) {
  UseMethod("BF", x)
}
