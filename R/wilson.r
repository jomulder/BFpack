#' Facial trustworthiness and criminal sentencing
#'
#' Data from a correlational study in which the correlation between ratings of
#' facial trustworthiness of inmates was correlated with whether they had
#' received the death penalty or not (wilson and Rule, 2015). These data were
#' simulated using the R-package \code{synthpop}, based on the characteristics
#' of the original data.
#'
#'
#'
#' \tabular{lll}{
#'    \strong{stim} \tab \code{integer} \tab Stimulus Number\cr
#'    \strong{sent} \tab \code{integer} \tab Sentence: 1 = Death, 0 = Life\cr
#'    \strong{race} \tab \code{integer} \tab Race: 1 = White, -1 = Black\cr
#'    \strong{glasses} \tab \code{integer} \tab Glasses: 1 = Yes, 0 = No\cr
#'    \strong{tattoos} \tab \code{integer} \tab Tattoos: 1 = Yes, 0 = No \cr
#'    \strong{ztrust} \tab \code{numeric} \tab Trustworthiness \cr
#'    \strong{trust_2nd} \tab \code{numeric} \tab Trustworthiness ratings with 2nd control group; Death targets are same as in primary analysis, Life targets are different.\cr
#'    \strong{afro} \tab \code{numeric} \tab raw Afrocentricity ratings.\cr
#'    \strong{zAfro} \tab \code{numeric} \tab Afrocentricity ratings normalized within target race. Analyses in paper were done with this variable.\cr
#'    \strong{attract} \tab \code{numeric} \tab Attractiveness\cr
#'    \strong{fWHR} \tab \code{numeric} \tab facial width-to-height \cr
#'    \strong{afWHR} \tab \code{numeric} \tab fWHR normalized within target race. Analyses in paper were done with this variable \cr
#'    \strong{maturity} \tab \code{numeric} \tab Maturity
#' }
#' @docType data
#' @keywords datasets
#' @name wilson
#' @usage data(wilson)
#' @references Wilson, J. P., & Rule, N. O. (2015). Facial Trustworthiness
#'   Predicts Extreme Criminal-Sentencing Outcomes. Psychological Science,
#'   26(8), 1325–1331. doi: 10.1177/0956797615590992
#' @format A data.frame with 742 rows and 13 columns.
NULL
