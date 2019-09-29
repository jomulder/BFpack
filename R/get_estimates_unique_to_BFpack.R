#' @method get_estimates hetcor
#' @import bain
#' @export
get_estimates.hetcor <- function(x, ...){
  cl <- match.call()
  cl[[1]] <- as.name("get_estimates")
  cl[["x"]] <- x$correlations
  P <- nrow(x$std.errors)
  out <- eval.parent(cl)
  stderr <- as.matrix(x$std.errors[lower.tri(diag(P))])
  errcov <- as.matrix(diag(stderr**2))
  out$Sigma <- list(errcov)
  out
}
