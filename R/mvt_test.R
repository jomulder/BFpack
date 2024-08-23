

#' @title Multivariate Student t test
#' @description First step to performs a Bayesian multivariate one sample Student t test using the
#' (adjusted) fractional Bayes factor using the \code{BF()} function.
#'
#'@details \code{X} must be a data matrix and \code{null}
#'must be a vector of the assumed null values of the variables.
#'
#'@param X a data matrix with the variables in the columns.
#'@param Y an optional data matrix with the variables in the columns.
#'@param null a vector of the null values of the variables.
#'@param paired a logical indicating whether you want a multivariate paired t-test.
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
mvt_test <- function(X, Y, null = NULL, paired = FALSE, ...){

  if(missing(Y)){

    X <- as.data.frame(as.matrix(X))
    p <- ncol(X)
    n <- nrow(X)

    if(is.null(null)){
      null <- rep(0,p)
    }

    if(length(null) != p){
      stop("'null' must be a vector of length equal to the number of variables.")
    }

    if(p == 1){
      out <- t_test(x=X,mu=null)
    }else{
      intercept <- rep(1,n)
      varnames <- colnames(X)
      formu <- as.formula(paste0("cbind(",paste0(varnames,collapse = ","),") ~ -1 + intercept"))
      out <- lm(formu,data=X)
      out$null <- null
      out$numpop <- 1
      out$paired <- paired
      class(out) <- "mvt_test"
      #
    }

  }else{

    if(paired==TRUE){

      p <- ncol(X)
      n <- nrow(X)

      if(is.null(null)){
        null <- rep(0,p)
      }

      if(length(null) != p){
        stop("'null' must be a vector of length equal to the number of variables.")
      }

      if(dim(X)[1]!=dim(Y)[1] || dim(X)[2]!=dim(Y)[2]){
        stop("X and Y must have same dimension in case of paired test")
      }

      if(p == 1){
        out <- t_test(x=X-Y,mu=null,paired=paired)
      }else{
        intercept <- rep(1,n)
        varnames <- colnames(X)
        formu <- as.formula(paste0("cbind(",paste0(varnames,collapse = ","),") ~ -1 + intercept"))
        df.XY <- as.data.frame(X-Y)
        colnames(df.XY) <- varnames
        out <- lm(formu,data=df.XY)
        out$null <- null
        out$numpop <- 1
        out$paired <- paired
        class(out) <- "mvt_test"
        #
      }

    }else{ # independent samples multivariate t test
      # equal covariances assumed

      p <- ncol(X)
      n1 <- nrow(X)
      n2 <- nrow(Y)
      N <- n1 + n2

      if(is.null(null)){
        null <- rep(0,p)
      }

      if(length(null) != p){
        stop("'null' must be a vector of length equal to the number of variables.")
      }

      if(dim(X)[2]!=dim(Y)[2]){
        stop("X and Y must have same number of columns")
      }

      if(p == 1){
        out <- t_test(x=X,y=Y,mu=null,paired=paired)
      }else{
        intercept <- rep(1,N)
        difference <- c(rep(1,n1),rep(0,n2))
        varnames <- colnames(X)
        formu <- as.formula(paste0("cbind(",paste0(varnames,collapse = ","),") ~ -1 + intercept + difference"))
        df.XY <- as.data.frame(rbind(X,Y))
        colnames(df.XY) <- varnames
        out <- lm(formu,data=df.XY)
        out$null <- null
        out$numpop <- 2
        out$paired <- paired
        class(out) <- "mvt_test"
        #
      }

    }

  }

  return(out)
}

#' @method get_estimates mvt_test
#' @export
get_estimates.mvt_test <- function(x, ...){
  class(x) <- "lm"
  x_est <- get_estimates(x)
  if(x$numpop==1){
    if(x$paired==FALSE){
      names(x_est$estimate) <- colnames(x_est$Sigma[[1]]) <-
        row.names(x_est$Sigma[[1]]) <- gsub("intercept_on_","",names(x_est$estimate))
    }else{
      names(x_est$estimate) <- colnames(x_est$Sigma[[1]]) <-
        row.names(x_est$Sigma[[1]]) <- gsub("intercept_on_","difference_",names(x_est$estimate))
    }
  }else{
    P <- ncol(x$coefficients)
    select_diff <- 2*(1:P)
    names(x_est$estimate) <- colnames(x_est$Sigma[[1]]) <-
      row.names(x_est$Sigma[[1]]) <- gsub("difference_on_","difference_",names(x_est$estimate))
    x_est$estimate <- x_est$estimate[select_diff]
    x_est$Sigma[[1]] <- x_est$Sigma[[1]][select_diff,select_diff]
  }
  x_est
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
                        BF.type = NULL,
                        ...) {

  if(is.null(BF.type)){
    BF.type <- "FBF"
    message("Note that the default Bayes factor has been changed to the fractional Bayes factor (FBF). To use
          the adjusted fractional Bayes factor (AFBF), set the 'BF.type' argument to 'AFBF'.")
  }

  if(x$numpop==1 & x$paired==FALSE){
    parameters <- "means"
  }else{
    parameters <- "differences"
  }

  if(x$numpop == 1){
    P <- length(x$coefficients)
    names1 <- colnames(x$coefficients)
    names2 <- paste0("intercept_on_",names1)
  }else{
    P <- ncol(x$coefficients)
    names1 <- colnames(x$coefficients)
    names2 <- paste0("difference_on_",names1)
  }
  #exploratory test of joint equality to null
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
  row.names(BF.explo$BFtu_confirmatory) <- row.names(BF.explo$PHP_confirmatory) <- parameters
  colnames(BF.explo$BFtu_confirmatory) <- c("BF0u","BFuu")
  colnames(BF.explo$PHP_confirmatory) <- c("Pr(=null)","Pr(not null)")

  if(!is.null(hypothesis)){
    variable.names <- colnames(x$coefficients)
    if(x$numpop==1 | x$paired==TRUE){
      add.name <- "intercept_on_"
      rem.name <- ""
      if(x$numpop==1 & x$paired==TRUE){
        rem.name <- "difference_"
      }
    }else{
      add.name <- "difference_on_"
      rem.name <- "difference_"
    }
    hypothesis.updated <- hypothesis
    for(varname in variable.names){
      hypothesis.updated <- gsub(paste0(rem.name,varname),paste0(add.name,varname),hypothesis.updated)
    }
  }else{
    hypothesis.updated <- NULL
  }

  BF.conf <- BF(x1,
                hypothesis=hypothesis.updated,
                prior.hyp.conf=prior.hyp.conf,
                log=log,
                complement=complement,
                BF.type=BF.type)
  if(!is.null(hypothesis)){
    BF.conf$hypotheses <- gsub(add.name,rem.name,BF.conf$hypotheses)
  }

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
    parameter=parameters,
    log = BF.conf$log,
    fraction_number_groupIDs = BF.conf$fraction_number_groupIDs,
    fraction_groupID_observations = BF.conf$fraction_groupID_observations,
    call=match.call())

  class(BFlm_out) <- "BF"

  return(BFlm_out)

}


