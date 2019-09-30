#BF method for rem.dyad class objects


#' @method BF rem.dyad
#' @export
BF.rem.dyad <- function(x,
                        hypothesis = NULL,
                        prior = NULL,
                        ...){

  #Extract summary statistics
  cl <- match.call()
  get_est <- get_estimates(x)
  cl[[1]] <- as.name("BF")
  cl[["x"]] <- get_est$estimate
  cl[["sigma"]] <- get_est$Sigma[[1]]
  cl[["n"]] <- x$m
  out <- eval.parent(cl)
  out$model <- x
  out$call <- match.call()
  out

}


