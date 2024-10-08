% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BF.mvt_test.R
\name{mvt_test}
\alias{mvt_test}
\title{Multivariate Student t test}
\usage{
mvt_test(X, Y, null = NULL, paired = FALSE, ...)
}
\arguments{
\item{X}{a data matrix with the variables in the columns.}

\item{Y}{an optional data matrix with the variables in the columns.}

\item{null}{a vector of the null values of the variables.}

\item{paired}{a logical indicating whether you want a multivariate paired t-test.}

\item{...}{further arguments to be passed to or from methods.}
}
\value{
An object that can be applied to the \code{BF()}.
}
\description{
First step to performs a Bayesian multivariate one sample Student t test using the
(adjusted) fractional Bayes factor using the \code{BF()} function.
}
\details{
\code{X} must be a data matrix and \code{null}
must be a vector of the assumed null values of the variables.
}
\examples{

mvt_fmri <- mvt_test(fmri[,1:2],null = c(0,0))
BF(mvt_fmri)

# the same test can be executed via the lm() function
intercept <- rep(1,nrow(fmri))
lm1 <- lm(cbind(Face,Vehicle) ~ -1 + intercept, data=fmri)
BF(lm1,hypothesis="intercept_on_Face=intercept_on_Vehicle=0")

}
\references{
Mulder, J. and Gu, X. (2023). Bayesian Testing of Scientific
Expectations under Multivariate Normal Linear Models. Multivariate Behavioral
Research, 57, 767-783. DOI: 10.1080/00273171.2021.1904809.
}
