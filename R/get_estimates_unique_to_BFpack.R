#' @method get_estimates hetcor
#' @import bain
#' @export
get_estimates.hetcor <- function(x, ...){
  cl <- match.call()
  cl[[1]] <- as.name("get_estimates.matrix")
  cl[["x"]] <- x$correlations
  eval.parent(cl)
}
