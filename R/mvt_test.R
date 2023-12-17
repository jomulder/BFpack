
#' @export
mvt_test <- function(Y, nullvalue = NULL, ...){

  Y <- as.matrix(Y)
  p <- ncol(Y)
  n <- nrow(Y)

  if(is.null(nullvalue)){
    nullvalue <- rep(0,p)
  }

  if(!missing(nullvalue[1]) && (length(nullvalue) != p || is.na(nullvalue[1])))
    stop("'nullvalue' must be a vector of length equal to the number of variables.")

  if(p == 1){
    out <- t_test(x=Y,mu=nullvalue)
  }else{
    intercept <- rep(1,n)
    varnames <- colnames(Y)
    out <- lm(as.formula(paste0("cbind(",paste0(varnames,collapse = ","),") ~ -1 + intercept")),data=Y)
    #
  }
  class(out) <- "mvt_test"

  return(out)
}

#' @method get_estimates glm
#' @export
get_estimates.mvt_test <- function(x, ...){
  out <- list()
  out$estimate <- coef(x)
  out$Sigma <- list(vcov(x))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "glm"
  out
}
