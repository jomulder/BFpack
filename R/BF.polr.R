#BF method for polr classes


#' @method BF polr
#' @export
BF.polr <- function(x,
                       hypothesis = NULL,
                       prior.hyp = NULL,
                       complement = TRUE,
                       BF.type = 2,
                       ...){

  #Extract summary statistics
  Args <- as.list(match.call()[-1])
  get_est <- get_estimates(x)
  Args$x <- get_est$estimate
  Args$Sigma <- get_est$Sigma[[1]]
  Args$n <- nrow(x$fitted.values)
  Args$hypothesis <- hypothesis
  Args$prior.hyp <- prior.hyp
  Args$complement <- complement
  Args$BF.type <- BF.type
  out <- do.call(BF, Args)
  out$model <- x
  out$call <- match.call()
  out

}


