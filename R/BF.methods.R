#' @title Bayes factors for Bayesian exploratory and confirmatory hypothesis
#' testing
#' @description The \code{BF} function can be used for hypothesis testing and
#'  model selection using the Bayes factor. By default exploratory hypothesis tests are
#' performed of whether each model parameter equals zero, is negative, or is
#'  positive. Confirmatory hypothesis tests can be executed by specifying hypotheses with
#' equality and/or order constraints on the parameters of interest. Depending on the
#' class of the fitted model different Bayes factors are used as described in the output.
#'
#' @param x An R object containing the outcome of a statistical analysis.
#' An R object containing the outcome of a statistical analysis. Currently, the
#' following objects can be processed: t_test(), bartlett_test(), lm(), aov(),
#' manova(), cor_test(), lmer() (only for testing random intercep variances),
#' glm(), coxph(), survreg(), polr(), zeroinfl(), rma(), ergm(), bergm(), or named
#' vector objects. In the case \code{x} is a named vector, the arguments \code{Sigma}
#' and \code{n} are also needed. See vignettes for elaborations.
#' @param hypothesis A character string containing the constrained (informative) hypotheses to
#' evaluate in a confirmatory test. The default is NULL, which will result in standard exploratory testing
#' under the model \code{x}.
#' @param prior.hyp.explo The prior probabilities of the hypotheses in the exploratory tests. Except for
#' objects of class \code{aov} (for (M)ANOVA, etc.), this argument should be a vector with three
#' elements reflecting the prior probability of a zero effect, a negative effect, and a positive effect, respectively.
#' For objects of class \code{aov}, the argument
#' should be a list where the first element should be a vector of length 3 specifying the prior probabilities of each parameter
#' being zero, negative, or positive, the second element should be a vector of length 2 specifying the prior probabilities
#' of a model where is no main effect for a factor and the full model, and the third element should be a vector of length 2
#' specifying the prior probabilities of a model where is no interaction effect (if present) for two factors and the full model.
#' The default (\code{NULL}) specifies equal prior probabilities for each hypothesis per exploratory test.
#' @param prior.hyp.conf The prior probabilities of the constrained hypotheses in the confirmatory test.
#' @param prior.hyp Deprecated. Please use the argument \code{prior.hyp.conf}.
#' @param complement a logical specifying whether the complement should be added
#' to the tested hypothesis under \code{hypothesis}.
#' @param log a logical specifying whether the Bayes factors should be computed on a log scale.
#' Default is \code{FALSE}.
#' @param BF.type An integer that specified the type of Bayes factor (or prior) that is used for the test.
#' Currently, this argument is only used for models of class 'lm' and 't_test',
#' where \code{BF.type=2} implies an adjusted fractional Bayes factor with a 'fractional prior mean' at the null value (Mulder, 2014),
#' and \code{BF.type=1} implies a regular fractional Bayes factor (based on O'Hagan (1995)) with a 'fractional prior mean' at the MLE.
#' @param iter Number of iterations that are used to compute the Monte Carlo estimates
#' (only used for certain hypotheses under multivariate models and when testing group variances).
#' @param Sigma An approximate posterior covariance matrix (e.g,. error covariance
#' matrix) of the parameters of interest. This argument is only required when \code{x}
#' is a named vector.
#' @param n The (effective) sample size that was used to acquire the estimates in the named vector
#' \code{x} and the error covariance matrix \code{Sigma}. This argument is only required when \code{x}
#' is a named vector.
#' @param ... Parameters passed to and from other functions.
#' @usage NULL
#' @return The output is an object of class \code{BF}. The object has elements:
#' \itemize{
#' \item \code{BFtu_exploratory}: The Bayes factors of the constrained hypotheses against
#' the unconstrained hypothesis in the exploratory test.
#' \item \code{BFtu_main} (only for \code{aov} objects with predictors of class \code{factor}):
#' The Bayes factors of a constrained model where all levels of a \code{factor} are assumed
#' to have the same effect on the outcome variable versus an unconstrained (full) model with
#' no constraints.
#' \item \code{BFtu_interaction} (only for \code{aov} objects with interaction effects with
#' predictors of class \code{factor}): The Bayes factors of a constrained model where the effect
#' of the dummy variables corresponding to an interaction effects are assumed to be zero versus
#' an unconstrained (full) model with no constraints.
#' \item \code{PHP_exploratory:} The posterior probabilities of the constrained hypotheses
#' in the exploratory test.
#' \item \code{PHP_main} (only for \code{aov} objects with predictors of class \code{factor}):
#' The posterior probabilities a constrained model where all levels of a \code{factor} are assumed
#' to have the same effect on the outcome variable versus an unconstrained (full) model with
#' no constraints.
#' \item \code{PHP_interaction} (only for \code{aov} objects with interaction effects with
#' predictors of class \code{factor}): The posterior probabilities of a constrained model where the
#' effect of the dummy variables corresponding to an interaction effects are assumed to be zero versus
#' an unconstrained (full) model with no constraints.
#' \item \code{BFtu_confirmatory}: The Bayes factors of the constrained hypotheses against
#' the unconstrained hypothesis in the confirmatory test using the \code{hypothesis}
#' argument.
#' \item \code{PHP_confirmatory}: The posterior probabilities of the constrained hypotheses
#' in the confirmatory test using the \code{hypothesis} argument.
#' \item \code{BFmatrix_confirmatory}: The evidence matrix which contains the Bayes factors
#' between all possible pairs of hypotheses in the confirmatory test.
#' \item \code{BFtable_confirmatory}: The \code{Specification table} (output when printing the
#' \code{summary} of a \code{BF} for a confirmatory test) which contains the different
#' elements of the extended Savage Dickey density ratio where
#' \itemize{
#' \item The first column `\code{complex=}' quantifies the relative complexity of the
#' equality constraints of a hypothesis (the prior density at the equality constraints in the
#' extended Savage Dickey density ratio).
#' \item The second column `\code{complex>}' quantifies the relative complexity of the
#' order constraints of a hypothesis (the prior probability of the order constraints in the extended
#' Savage Dickey density ratio).
#' \item The third column `\code{fit=}' quantifies the relative fit of the equality
#' constraints of a hypothesis (the posterior density at the equality constraints in the extended
#' Savage Dickey density ratio).
#' \item The fourth column `\code{fit>}' quantifies the relative fit of the order
#' constraints of a hypothesis (the posterior probability of the order constraints in the extended
#' Savage Dickey density ratio)
#' \item The fifth column `\code{BF=}' contains the Bayes factor of the equality constraints
#' against the unconstrained hypothesis.
#' \item The sixth column `\code{BF>}' contains the Bayes factor of the order constraints
#' against the unconstrained hypothesis.
#' \item The seventh column `\code{BF}' contains the Bayes factor of the constrained hypothesis
#' against the unconstrained hypothesis.
#' \item The eighth column `\code{PHP}' contains the posterior probabilities of the hypotheses.
#' }
#' \item \code{prior.hyp.explo}: The prior probabilities of the constrained hypotheses in the exploratory tests.
#' \item \code{prior.hyp.conf}: The prior probabilities of the constrained hypotheses in the confirmatory test.
#' \item \code{hypotheses}: The tested constrained hypotheses in a confirmatory test.
#' \item \code{estimates}: The unconstrained estimates.
#' \item \code{model}: The input model \code{x}.
#' \item \code{bayesfactor}: The type of Bayes factor that is used for this model.
#' \item \code{parameter}: The type of parameter that is tested.
#' \item \code{log}: \code{logical} whether the Bayes factors were reported on a log scale.
#' \item \code{fraction_number_groupIDs} (only for objects of class \code{lm}): The number of
#' 'group identifiers' that were identified based on the number of unique combinations of \code{level}s
#' of predictor variables of class \code{factor} in the data. These group identifiers are used to automatically
#' specify the minimal fractions that are used to compute (adjusted) fractional Bayes factors.
#' \item \code{fraction_groupID_observations} (only for objects of class \code{lm}): A vector that
#' specifies to which 'group identifier' an observation belongs. The group identifiers are constructed
#' based on the unique combination of the \code{levels} based on the predictor variables of class \code{factor}
#' of the observations.
#' \item \code{call}: The call of the \code{BF} function.
#' }
#' @details The function requires a fitted modeling object. Current analyses
#' that are supported: \code{\link[bain]{t_test}},
#' \code{\link[BFpack]{bartlett_test}},
#' \code{\link[stats]{aov}}, \code{\link[stats]{manova}},
#' \code{\link[stats]{lm}}, \code{mlm},
#' \code{\link[stats]{glm}}, \code{\link[polycor]{hetcor}},
#' \code{\link[lme4]{lmer}}, \code{\link[survival]{coxph}},
#' \code{\link[survival]{survreg}}, \code{\link[ergm]{ergm}},
#' \code{\link[Bergm]{bergm}},
#' \code{\link[pscl]{zeroinfl}}, \code{\link[metafor]{rma}} and \code{\link[MASS]{polr}}.
#'
#' For testing parameters from the results of t_test(), lm(), aov(),
#' manova(), and bartlett_test(), hypothesis testing is done using
#' adjusted fractional Bayes factors are computed (using minimal fractions).
#' For testing measures of association (e.g., correlations) via \code{cor_test()},
#' Bayes factors are computed using joint uniform priors under the correlation
#' matrices. For testing intraclass correlations (random intercept variances) via
#' \code{lmer()}, Bayes factors are computed using uniform priors for the intraclass
#' correlations. For all other tests,  approximate adjusted fractional Bayes factors
#' (with minimal fractions) are computed using Gaussian approximations, similar as
#' a classical Wald test.
#'
#' @references Mulder, J., D.R. Williams, Gu, X., A. Tomarken,
#' F. Böing-Messing, J.A.O.C. Olsson-Collentine, Marlyne Meyerink, J. Menke,
#' J.-P. Fox, Y. Rosseel, E.J. Wagenmakers, H. Hoijtink., and van Lissa, C.
#' (2021). BFpack: Flexible Bayes Factor Testing of Scientific Theories
#' in R. Journal of Statistical Software. <https://doi.org/10.18637/jss.v100.i18>
#'
#' @examples
#' \donttest{
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
#' # EXAMPLE 3. linear regression
#' lm1 <- lm(mpg ~ cyl + hp + wt, data = mtcars)
#' BF(lm1, hypothesis = "wt < cyl < hp = 0")
#'
#' # EXAMPLE 4. Logistic regression
#' fit <- glm(sent ~ ztrust + zfWHR + zAfro + glasses + attract + maturity +
#'    tattoos, family = binomial(), data = wilson)
#' BF1 <- BF(fit, hypothesis = "ztrust > zfWHR > 0;
#'                              ztrust > 0 & zfWHR = 0")
#' summary(BF1)
#'
#' # EXAMPLE 5. Correlation analysis
#' #set.seed(123)
#' #cor1 <- cor_test(memory[1:20,1:3])
#' #BF1 <- BF(cor1)
#' #summary(BF1)
#' #BF2 <- BF(cor1, hypothesis = "Wmn_with_Im > Wmn_with_Del > 0;
#' #                             Wmn_with_Im = Wmn_with_Del = 0")
#' #summary(BF2)
#'
#' # EXAMPLE 6. Bayes factor testing on a named vector
#' # A Poisson regression model is used to illustrate the computation
#' # of Bayes factors with a named vector as input
#' poisson1 <- glm(formula = breaks ~ wool + tension,
#'   data = datasets::warpbreaks, family = poisson)
#' # extract estimates, error covariance matrix, and sample size:
#' estimates <- poisson1$coefficients
#' covmatrix <- vcov(poisson1)
#' samplesize <- nobs(poisson1)
#' # compute Bayes factors on equal/order constrained hypotheses on coefficients
#' BF1 <- BF(estimates, Sigma = covmatrix, n = samplesize, hypothesis =
#' "woolB > tensionM > tensionH; woolB = tensionM = tensionH")
#' summary(BF1)
#' }
#' @rdname BF
#' @export
#' @useDynLib BFpack, .registration = TRUE
#'
BF <- function(x, hypothesis, prior.hyp.explo, prior.hyp.conf, prior.hyp, complement, log, ...) {
  UseMethod("BF", x)
}
