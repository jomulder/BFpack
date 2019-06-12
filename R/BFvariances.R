#' @title Bartlett Test of Homogeneity of Variances
#' @description Performs Bartlett's test of the null that the variances in each of the groups (samples) are the same.
#' @param x a numeric vector of data values, or a list of numeric data vectors representing the respective samples, or fitted linear model objects (inheriting from class "\code{lm}").
#' @param g a vector or factor object giving the group for the corresponding elements of \code{x}. Ignored if \code{x} is a list.
#' @param formula a formula of the form \code{lhs ~ rhs} where lhs gives the data values and rhs the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see \code{\link{model.frame}}) containing the variables in the formula \code{formula}. By default the variables are taken from environment(\code{formula}).
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when the data contain \code{NA}s. Defaults to getOption("\code{na.action}").
#' @section BF var_test:
#' In order to allow users to enjoy the functionality of BFpack with the familiar
#' stats-function bartlett.test, we have had to make minor changes to the function
#' bartlett.test.default. All rights to, and credit for, the function bartlett.test.default
#' belong to the R Core Team, as indicated in the original license below.
#' We make no claims to copyright and incur no liability with regard to the
#' changes implemented in var_test.
#'
#' This the original copyright notice by the R core team:
#' File src/library/stats/R/t_test.R
#' Part of the R package, https://www.R-project.org
#'
#' Copyright (C) 1995-2015 The R Core Team
#'
#' This program is free software; you can redistribute it and/or modify
#' it under the terms of the GNU General Public License as published by
#' the Free Software Foundation; either version 2 of the License, or
#' (at your option) any later version.
#'
#' This program is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU General Public License for more details.
#'
#' A copy of the GNU General Public License is available at
#' https://www.R-project.org/Licenses/
#'
#' @aliases var_test var_test.default var_test.formula
#'
#' @return A list with class \code{"vtest"} containing the following
#' components:
#' \item{statistic}{Bartlett's K-squared test statistic.}
#' \item{parameter}{the degrees of freedom of the approximate chi-squared distribution of the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{method}{the character string "\code{Bartlett test of homogeneity of variances}".}
#' \item{data.name}{a character string giving the name(s) of
#' the data.}
#' \item{v}{a vector containing the sample variances for each factor level.}
#' \item{n}{a vector containing the sample sizes for each factor level.}
#' @details If \code{x} is a list, its elements are taken as the samples or fitted linear models to be compared for homogeneity of variances. In this case, the elements must either all be numeric data vectors or fitted linear model objects, \code{g} is ignored, and one can simply use \code{bartlett.test(x)} to perform the test. If the samples are not yet contained in a list, use \code{var_test(list(x, ...))}.
#'
#' Otherwise, x must be a numeric data vector, and g must be a vector or factor object of the same length as x giving the group for the corresponding elements of x.
#' @references Bartlett, M. S. (1937). Properties of sufficiency and statistical tests. Proceedings of the Royal Society of London Series A 160, 268–282. doi: 10.1098/rspa.1937.0109.
#' @seealso \code{\link{var.test}} for the special case of comparing variances in two samples from normal distributions; \code{\link{fligner.test} for a rank-based (nonparametric) k-sample test for homogeneity of variances; \code{\link{ansari.test} and \code{\link{mood.test} for two rank based two-sample tests for difference in scale.
#' @examples
#' require(graphics)
#'
#' plot(count ~ spray, data = InsectSprays)
#' var_test(InsectSprays$count, InsectSprays$spray)
#' var_test(count ~ spray, data = InsectSprays)
#'
#' @rdname var_test
#' @export
var_test <- function (x, ...) UseMethod("var_test")

#' @method var_test default
#' @rdname var_test
#' @export
var_test.default <- function (x, g, ...)
{
  LM <- FALSE
  if (is.list(x)) {
    if (length(x) < 2L)
      stop("'x' must be a list with at least 2 elements")
    DNAME <- deparse(substitute(x))
    if (all(sapply(x, function(obj) inherits(obj, "lm"))))
      LM <- TRUE
    else x <- lapply(x, function(x) x <- x[is.finite(x)])
    k <- length(x)
  }
  else {
    if (length(x) != length(g))
      stop("'x' and 'g' must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    OK <- complete.cases(x, g)
    x <- x[OK]
    g <- factor(g[OK])
    k <- nlevels(g)
    if (k < 2)
      stop("all observations are in the same group")
    x <- split(x, g)
  }
  if (LM) {
    n <- sapply(x, function(obj) obj$df.resid)
    v <- sapply(x, function(obj) sum(obj$residuals^2))/n
  }
  else {
    n <- sapply(x, "length") - 1
    if (any(n <= 0))
      stop("there must be at least 2 observations in each group")
    v <- sapply(x, "var")
  }
  n.total <- sum(n)
  v.total <- sum(n * v)/n.total
  STATISTIC <- ((n.total * log(v.total) - sum(n * log(v)))/(1 +
                                                              (sum(1/n) - 1/n.total)/(3 * (k - 1))))
  PARAMETER <- k - 1
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  names(STATISTIC) <- "Bartlett's K-squared"
  names(PARAMETER) <- "df"
  RVAL <- list(statistic = STATISTIC, parameter = PARAMETER,
               p.value = PVAL, data.name = DNAME, method = "Bartlett test of homogeneity of variances",
               n = n, v = v)
  class(RVAL) <- "vtest"
  return(RVAL)
}

#' @method var_test formula
#' @rdname var_test
#' @export
var_test.formula <- function (formula, data, subset, na.action, ...)
{
  if (missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  mf <- eval(m, parent.frame())
  if (length(mf) != 2L)
    stop("'formula' should be of the form response ~ group")
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  y <- do.call("var_test", as.list(mf))
  y$data.name <- DNAME
  y
}


#' @method BF vtest
#' @rdname BF
#' @export
BF.vtest <- function(x,
                     hypothesis = NULL,
                     prior = NULL,
                     parameter = NULL,
                     ...){
  vars <- x$v
  nvec <- x$n
  var_names <- names(x$v)

  #####
  #
  # FLORIANS FUNCTIONS
  #
  #####

}



