#BF method for survreg classes


#' @method BF survreg
#' @export
BF.survreg <- function(x,
                   hypothesis = NULL,
                   prior.hyp = NULL,
                   complement = TRUE,
                   log = FALSE,
                   ...){

  #Extract summary statistics
  get_est <- get_estimates(x)
  Args <- as.list(match.call()[-1])
  get_est <- get_estimates(x)
  Args$x <- get_est$estimate
  Args$Sigma <- get_est$Sigma[[1]]
  Args$n <- length(x$y)
  Args$hypothesis <- hypothesis
  Args$prior.hyp <- prior.hyp
  Args$complement <- complement
  Args$log <- log
  out <- do.call(BF, Args)
  out$model <- x
  out$call <- match.call()
  out
}


