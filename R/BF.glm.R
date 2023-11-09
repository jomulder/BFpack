#BF method for glm classes


#' @method BF glm
#' @export
BF.glm <- function(x,
                   hypothesis = NULL,
                   prior.hyp = NULL,
                   complement = TRUE,
                   log = FALSE,
                   ...){

  Args <- as.list(match.call()[-1])
  get_est <- get_estimates(x)
  Args$x <- get_est$estimate
  Args$Sigma <- get_est$Sigma[[1]]
  Args$n <- nobs(x)
  Args$hypothesis <- hypothesis
  Args$prior.hyp <- prior.hyp
  Args$complement <- complement
  Args$log <- log
  out <- do.call(BF, Args)
  out$model <- x
  out$call <- match.call()
  out
}


