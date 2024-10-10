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


#' @method get_estimates lm
#' @export
get_estimates.lm <- function(x, ...){
  out <- list()
  P <- ncol(x$coefficients)
  K <- nrow(x$coefficients)
  N <- nrow(x$residuals)
  if(!is.matrix(x$coefficients)){
    out$estimate <- coef(x)
    out$Sigma <- list(vcov(x))
    class(out) <- "model_estimates"
    attr(out, "analysisType") <- "lm"
    out
  }else{
    names_coef1 <- row.names(x$coefficients)
    names_coef2 <- colnames(x$coefficients)
    names_coef <- unlist(lapply(1:P,function(p){
      lapply(1:K,function(k){
        paste0(names_coef1[k],"_on_",names_coef2[p])
      })
    }))
    # estimates of regression coefficients
    estimatesBeta <- c(x$coefficients)
    names(estimatesBeta) <- names_coef
    Xmat <- model.matrix(x)
    Ymat <- model.matrix(x)%*%x$coefficients + x$residuals
    SigmaEst <- t(x$residuals)%*%x$residuals/N
    covmBeta <- kronecker(SigmaEst,solve(t(Xmat)%*%Xmat))
    row.names(covmBeta) <- colnames(covmBeta) <- names_coef

    out$estimate <- estimatesBeta
    out$Sigma <- list(covmBeta)
    class(out) <- "model_estimates"
    attr(out, "analysisType") <- "mlm"
    out
  }
}


#' @method get_estimates cor_test
#' @export
get_estimates.cor_test <- function(x, ...){
  out <- list()
  out$estimate <- x$meanF
  out$Sigma <- list(x$covmF)
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "corr_htest"
  out
}


#' @method get_estimates t_test
#' @export
get_estimates.t_test <- function(x, ...){
  out <- list()
  if(length(x$estimate)>1){
    difference <- x$estimate[1] - x$estimate[2]
    names(difference) <- "difference"
    out$estimate <- difference
    out$Sigma <- list((x$stderr)**2)
  }else if(names(x$estimate) == "mean difference"){
    difference <- x$estimate
    names(difference) <- "difference"
    out$estimate <- difference
    out$Sigma <- list((x$stderr)**2)
  }else{
    populationmean <- x$estimate
    names(populationmean) <- "mu"
    out$estimate <- populationmean
    out$Sigma <- list((x$stderr)**2)
  }
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "t_test"
  out
}




