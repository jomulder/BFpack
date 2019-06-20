#' @title BF test for a univariate normal linear regression type model
#' @description FUNCTION_DESCRIPTION
#' @param x An R object containing the outcome of a statistical analysis.
#' @param hypothesis A character string containing the informative hypotheses to
#' evaluate. Defaults to NULL, which will
#' result in an exploratory analysis.
#' @param prior Prior, Default: 'default'
#' @param parameter Joris, wat is dit?
#' @param ... Parameters passed to and from other functions.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname BF
#' @export
#' @useDynLib BFpack, .registration = TRUE
#'
# The functions below can be used for Bayesian hypothesis testing under univariate
# and multivariate linear regression models. Also grouping variables (factors) are
# properly dealt with using group specific fractions for the default priors.
BF <- function(x, hypothesis, prior, parameter, ...) {
  UseMethod("BF", x)
}
