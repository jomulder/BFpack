#BF method for coeftest class objects

#' @method BF coeftest
#' @export
BF.coeftest <- function(x,
                     hypothesis = NULL,
                     prior.hyp.explo = NULL,
                     prior.hyp.conf = NULL,
                     prior.hyp = NULL,
                     complement = TRUE,
                     log = FALSE,
                     ...){

  logIN <- log

  Sigma <- diag(x[, 2L]^2)
  n <- attr(x, "nobs")

  if(is.null(n)) stop("'BF.coeftest' only works if 'nobs.coeftest' gives the number of observations.")
  if(!is.null(hypothesis)) warning("constrained hypothesis testing is not supported for objects of class 'coeftest'")
  if(!is.null(prior.hyp)) warning("prior specification via 'prior.hyp' is not supported for objects of class 'coeftest'")
  #if(!exploratory) stop("only exploratory hypothesis testing is supported for objects of class 'coeftest'")

  out <- BF.default(x[, 1L], Sigma = Sigma, n = n, log = logIN, ...)
  out$model <- x
  out$call <- match.call()
  out

}




