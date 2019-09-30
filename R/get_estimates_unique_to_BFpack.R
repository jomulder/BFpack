#' @method get_estimates hetcor
#' @import bain
#' @export
get_estimates.hetcor <- function(x, ...){
  cl <- match.call()
  cl[[1]] <- as.name("get_estimates")
  cl[["x"]] <- x$correlations
  P <- nrow(x$std.errors)
  out <- eval.parent(cl)
  retain <- matrix(1:length(out$estimate), nrow = nrow(x$std.errors))
  out$estimate <- out$estimate[retain[lower.tri(retain)]]
  errcov <- x$std.errors**2
  errcov <- errcov[lower.tri(retain)]
  if(length(errcov) == 1){
    out$Sigma <- list(matrix(errcov))
  } else {
    out$Sigma <- list(diag(errcov))
  }
  out
}
