#BF method for polr classes


#' @method BF polr
#' @export
BF.polr <- function(x,
                       hypothesis = NULL,
                       prior.hyp.explo = NULL,
                       prior.hyp.conf = NULL,
                       prior.hyp = NULL,
                       complement = TRUE,
                       log = FALSE,
                    cov.prob = .95,
                       ...){

  #Extract summary statistics
  Args <- as.list(match.call()[-1])
  get_est <- get_estimates(x)
  Args$x <- get_est$estimate
  Args$Sigma <- get_est$Sigma[[1]]
  Args$n <- nrow(x$fitted.values)
  Args$hypothesis <- hypothesis
  Args$prior.hyp <- prior.hyp
  Args$prior.hyp.explo <- prior.hyp.explo
  Args$prior.hyp.conf <- prior.hyp.conf
  Args$complement <- complement
  Args$log <- log
  Args$cov.prob <- cov.prob
  out <- do.call(BF, Args)
  out$model <- x
  out$call <- match.call()
  out

}


