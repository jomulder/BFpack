

#' @title Multivariate Student t test
#' @description First step to performs a Bayesian multivariate one sample Student t test using the
#' (adjusted) fractional Bayes factor using the \code{BF()} function.
#'
#'@details \code{Y} must be a data matrix and \code{null}
#'must be a vector of the assumed null values of the variables.
#'
#'@param Y a data matrix with different variables in the columns.
#'@param null a vector of the null values of the variables.
#'@param ... further arguments to be passed to or from methods.
#'
#'@return An object that can be applied to the \code{BF()}.
#'
#'@references Mulder, J. and Gu, X. (2023). Bayesian Testing of Scientific
#'Expectations under Multivariate Normal Linear Models. Multivariate Behavioral
#'Research, 57, 767-783. DOI: 10.1080/00273171.2021.1904809.
#'
#'@examples
#'
#'mvt_fmri <- mvt_test(fmri[,1:2],null = c(0,0))
#'BF(mvt_fmri)
#'
#'# the same test can be executed via the lm() function
#'intercept <- rep(1,nrow(fmri))
#'lm1 <- lm(cbind(Face,Vehicle) ~ -1 + intercept, data=fmri)
#'BF(lm1,hypothesis="intercept_on_Face=intercept_on_Vehicle=0")
#'
#' @rdname mvt_test
#' @export
mvt_test <- function(Y, null = NULL, ...){

  Y <- as.data.frame(as.matrix(Y))
  p <- ncol(Y)
  n <- nrow(Y)

  if(is.null(null)){
    null <- rep(0,p)
  }

  if(length(null) != p){
    stop("'null' must be a vector of length equal to the number of variables.")
  }

  if(p == 1){
    out <- t_test(x=Y,mu=null)
  }else{
    intercept <- rep(1,n)
    varnames <- colnames(Y)
    formu <- as.formula(paste0("cbind(",paste0(varnames,collapse = ","),") ~ -1 + intercept"))
    out <- lm(formu,data=Y)
    out$null <- null
    class(out) <- "mvt_test"
    #
  }

  return(out)
}

#' @method get_estimates mvt_test
#' @export
get_estimates.mvt_test <- function(x, ...){
  class(x) <- "lm"
  get_estimates(x)
}

#' @method BF mvt_test
#' @export
BF.mvt_test <- function(x,
                        hypothesis = NULL,
                        prior.hyp.explo = NULL,
                        prior.hyp.conf = NULL,
                        prior.hyp = NULL,
                        complement = TRUE,
                        log = FALSE,
                        BF.type = 2,
                        ...) {

  P <- length(x$coefficients)
  names1 <- colnames(x$coefficients)
  names2 <- paste0("intercept_on_",names1)

  hypothesis.explo <- paste0(unlist(lapply(1:P,function(p){
    paste0(names2[p],"=",x$null[p])
  })),collapse = " & ")
  x1 <- x
  class(x1) <- "lm"
  BF.explo <- BF(x1,
                 hypothesis=hypothesis.explo,
                 prior.hyp.conf=prior.hyp.explo,
                 log=log,
                 BF.type=BF.type)
  BF.explo$BFtu_confirmatory <- t(as.matrix(BF.explo$BFtu_confirmatory))
  BF.explo$PHP_confirmatory <- t(as.matrix(BF.explo$PHP_confirmatory))
  row.names(BF.explo$BFtu_confirmatory) <- row.names(BF.explo$PHP_confirmatory) <- "means"
  colnames(BF.explo$BFtu_confirmatory) <- c("BF0u","BFuu")
  colnames(BF.explo$PHP_confirmatory) <- c("Pr(=null)","Pr(not null)")

  BF.conf <- BF(x1,
                hypothesis=hypothesis,
                prior.hyp.conf=prior.hyp.conf,
                log=log,
                complement=complement,
                BF.type=BF.type)

  BFlm_out <- list(
    BFtu_exploratory=BF.explo$BFtu_confirmatory,
    PHP_exploratory=BF.explo$PHP_confirmatory,
    BFtu_confirmatory=BF.conf$BFtu_confirmatory,
    PHP_confirmatory=BF.conf$PHP_confirmatory,
    BFmatrix_confirmatory=BF.conf$BFmatrix_confirmatory,
    BFtable_confirmatory=BF.conf$BFtable_confirmatory,
    prior.hyp.explo=BF.explo$prior.hyp.conf,
    prior.hyp.conf=BF.conf$prior.hyp.conf,
    hypotheses=BF.conf$hypotheses,
    estimates=BF.conf$estimates,
    model=x1,
    bayesfactor=BF.conf$bayesfactor,
    parameter="means",
    log = BF.conf$log,
    fraction_number_groupIDs = BF.conf$fraction_number_groupIDs,
    fraction_groupID_observations = BF.conf$fraction_groupID_observations,
    call=match.call())

  class(BFlm_out) <- "BF"

  return(BFlm_out)

}


