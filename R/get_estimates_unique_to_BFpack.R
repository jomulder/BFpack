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
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "hetcor"
  out
}


#' @method get_estimates coxph
#' @export
get_estimates.coxph <- function(x, ...){
  out <- list()
  out$estimate <- coef(x)
  out$Sigma <- list(vcov(x))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "cophx"
  out
}


#' @method get_estimates glm
#' @export
get_estimates.glm <- function(x, ...){
  out <- list()
  out$estimate <- coef(x)
  out$Sigma <- list(vcov(x))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "glm"
  out
}

#' @method get_estimates polr
#' @export
get_estimates.polr <- function(x, ...){
  out <- list()
  out$estimate <- c(coef(x),x$zeta)
  out$Sigma <- list(vcov(x))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "polr"
  out
}


#' @method get_estimates bartlett_htest
#' @export
get_estimates.bartlett_htest <- function(x, ...){
  out <- list()
  out$estimate <- x$vars
  out$Sigma <- NULL
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "bartlett_htest"
  out
}


#' @method get_estimates survreg
#' @export
get_estimates.survreg <- function(x, ...){
  out <- list()
  out$estimate <- x$coefficients
  out$Sigma <- list(x$var)
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "survreg"
  out
}


#' @method get_estimates zeroinfl
#' @export
get_estimates.zeroinfl <- function(x, ...){
  out <- list()
  out$estimate <- c(coef(x),x$zeta)
  out$Sigma <- list(vcov(x))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "zeroinfl"
  out
}





