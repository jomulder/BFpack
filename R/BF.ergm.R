#BF method for ergm class

#' @importFrom sandwich sandwich
#' @importFrom ergm ergmMPLE
#' @method BF ergm
#' @export
BF.ergm <- function(x,
                    hypothesis = NULL,
                    prior.hyp = NULL,
                    complement = TRUE,
                    ...){

  # extract coefficients
  estimate <- coef(x)
  K1 <- length(estimate)
  # get design matrix of pseudo likelihood to construct prior covariance matrix
  x_MPLE <- ergm::ergmMPLE(formula=x$formula,output="dyadlist")
  design.X <- x_MPLE$predictor[,2+1:K1]
  which.edges <- which(colnames(design.X)=="edges")
  if(length(which.edges)==0){ #no intercept 'edges'
    Xdelta <- as.matrix(design.X)
    priorcov <- solve(t(Xdelta)%*%Xdelta) * nrow(Xdelta)
  }else{
    Xdelta <- as.matrix(design.X[,-which.edges])
    priorcov.Xdelta <- solve(t(Xdelta)%*%Xdelta) * nrow(Xdelta)
    priorcov <- matrix(0,ncol=K1,nrow=K1)
    priorcov[which.edges,which.edges] <- 100000 #flat prior for the 'edges' parameter
    if(which.edges==1){
      priorcov[2:K1,2:K1] <- priorcov.Xdelta
    }else{
      if(which.edges==K1){
        priorcov[1:(K1-1),1:(K1-1)] <- priorcov.Xdelta
      }else{
        priorcov[1:(which.edges-1),1:(which.edges-1)] <- priorcov.Xdelta[1:(which.edges-1),1:(which.edges-1)]
        priorcov[(which.edges+1):K1,(which.edges+1):K1] <- priorcov.Xdelta[which.edges:(K1-1),which.edges:(K1-1)]
        priorcov[1:(which.edges-1),(which.edges+1):K1] <- priorcov.Xdelta[1:(which.edges-1),which.edges:(K1-1)]
        priorcov[(which.edges+1):K1,1:(which.edges-1)] <- t(priorcov[1:(which.edges-1),(which.edges+1):K1])
      }
    }
  }
  Bergm.out <- Bergm::bergm(x$formula,prior.mean=rep(0,K1),prior.sigma=priorcov,...)
  #get robust estimates for the Gaussian mean and covariance matrix
  post.mean <- apply(Bergm.out$Theta,2,median)
  names(post.mean) <- names(estimate)
  #get robust estimate of posterior covariance matrix
  mlm1 <- lm(Bergm.out$Theta ~ 1)
  post.Sigma <- sandwich(mlm1) * nrow(Bergm.out$Theta)

  # use Savage-Dickey approximation of the BF
  if(length(which.edges)==0){
    prior.mean = rep(0,K1)
    names(prior.mean) <- names(estimate)
    row.names(priorcov) <- colnames(priorcov) <- names(estimate)
    BFergm_out <- Savage.Dickey.Gaussian(prior.mean = prior.mean,
                     prior.sigma = priorcov,
                     post.mean = post.mean,
                     post.sigma = post.Sigma,
                     hypothesis = hypothesis,
                     prior.hyp = prior.hyp,
                     complement = complement)
  }else{
    prior.mean = rep(0,K1)
    names(prior.mean) <- names(estimate)
    row.names(priorcov) <- colnames(priorcov) <- names(estimate)
    BFergm_out <- Savage.Dickey.Gaussian(prior.mean = prior.mean[-which.edges],
                                         prior.sigma = priorcov[-which.edges,-which.edges],
                                         post.mean = post.mean[-which.edges],
                                         post.sigma = post.Sigma[-which.edges,-which.edges],
                                         hypothesis = hypothesis,
                                         prior.hyp = prior.hyp,
                                         complement = complement)
  }

  BFergm_out$model <- x
  BFergm_out$call <- match.call()
  BFergm_out$bayesfactor <- "Bayes factors based on unit-information priors and Gaussian approximations"
  BFergm_out$parameter <- "ERGM coefficients"
  BFergm_out$model_update <- Bergm.out
  BFergm_out$prior.parameters <- list(prior.mean=prior.mean,prior.cov=priorcov)

  return(BFergm_out)
}

#' @method get_estimates ergm
#' @export
get_estimates.ergm <- function(x, ...){
  estimate <- coef(x)
  K1 <- length(estimate)
  x_MPLE <- ergm::ergmMPLE(formula=x$formula,output="dyadlist")
  design.X <- x_MPLE$predictor[,2+1:K1]
  which.edges <- which(colnames(design.X)=="edges")
  out <- list()
  if(length(which.edges)==0){ #no intercept 'edges'
    out$estimate <- coef(x)
    out$Sigma <- list(vcov(x))
  }else{
    out$estimate <- coef(x)[-which.edges]
    out$Sigma <- list(as.matrix(vcov(x)[-which.edges,-which.edges]))
  }
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "ergm"
  out
}


#' @method BF bergm
#' @export
BF.bergm <- function(x,
                    hypothesis = NULL,
                    prior.hyp = NULL,
                    complement = TRUE,
                    ...){

  # first check if effect names in hypothesis argument correspond with names in x
  coef_names_hyp <- names(get_estimates(x)$estimate)
  if(!is.null(hypothesis)){
    test0 <- parse_hypothesis(coef_names_hyp,hypothesis)
  }

  # extract coefficients
  estimate <- apply(x$Theta,2,median)
  K1 <- length(estimate)
  coef_names <- paste0("theta",1:K1)
  # get design matrix of pseudo likelihood to construct prior covariance matrix
  x_MPLE <- ergm::ergmMPLE(formula=x$formula,output="dyadlist")
  design.X <- x_MPLE$predictor[,2+1:K1]
  which.edges <- which(colnames(design.X)=="edges")
  if(length(which.edges)==0){ #no intercept 'edges'
    Xdelta <- as.matrix(design.X)
    priorcov <- solve(t(Xdelta)%*%Xdelta) * nrow(Xdelta)
  }else{
    Xdelta <- as.matrix(design.X[,-which.edges])
    priorcov.Xdelta <- solve(t(Xdelta)%*%Xdelta) * nrow(Xdelta)
    priorcov <- matrix(0,ncol=K1,nrow=K1)
    priorcov[which.edges,which.edges] <- 100000 #flat prior for the 'edges' parameter
    if(which.edges==1){
      priorcov[2:K1,2:K1] <- priorcov.Xdelta
    }else{
      if(which.edges==K1){
        priorcov[1:(K1-1),1:(K1-1)] <- priorcov.Xdelta
      }else{
        priorcov[1:(which.edges-1),1:(which.edges-1)] <- priorcov.Xdelta[1:(which.edges-1),1:(which.edges-1)]
        priorcov[(which.edges+1):K1,(which.edges+1):K1] <- priorcov.Xdelta[which.edges:(K1-1),which.edges:(K1-1)]
        priorcov[1:(which.edges-1),(which.edges+1):K1] <- priorcov.Xdelta[1:(which.edges-1),which.edges:(K1-1)]
        priorcov[(which.edges+1):K1,1:(which.edges-1)] <- t(priorcov[1:(which.edges-1),(which.edges+1):K1])
      }
    }
  }
  Bergm.out <- Bergm::bergm(x$formula,prior.mean=rep(0,K1),prior.sigma=priorcov,...)
  #get robust estimates for the Gaussian mean and covariance matrix
  post.mean <- apply(Bergm.out$Theta,2,median)
  names(post.mean) <- paste0("theta",1:K1)
  #get robust estimate of posterior covariance matrix
  mlm1 <- lm(Bergm.out$Theta ~ 1)
  post.Sigma <- sandwich(mlm1) * nrow(Bergm.out$Theta)

  # use Savage-Dickey approximation of the BF
  if(length(which.edges)==0){
    prior.mean = rep(0,K1)
    names(prior.mean) <- coef_names
    row.names(priorcov) <- colnames(priorcov) <- coef_names
    BFergm_out <- Savage.Dickey.Gaussian(prior.mean = prior.mean,
                                         prior.sigma = priorcov,
                                         post.mean = post.mean,
                                         post.sigma = post.Sigma,
                                         hypothesis = hypothesis,
                                         prior.hyp = prior.hyp,
                                         complement = complement)
  }else{
    prior.mean = rep(0,K1)
    names(prior.mean) <- coef_names
    row.names(priorcov) <- colnames(priorcov) <- coef_names
    BFergm_out <- Savage.Dickey.Gaussian(prior.mean = prior.mean[-which.edges],
                                         prior.sigma = priorcov[-which.edges,-which.edges],
                                         post.mean = post.mean[-which.edges],
                                         post.sigma = post.Sigma[-which.edges,-which.edges],
                                         hypothesis = hypothesis,
                                         prior.hyp = prior.hyp,
                                         complement = complement)
  }

  BFergm_out$model <- x
  BFergm_out$call <- match.call()
  BFergm_out$bayesfactor <- "Bayes factors based on unit-information priors and Gaussian approximations"
  BFergm_out$parameter <- "ERGM coefficients"
  BFergm_out$model_update <- Bergm.out
  BFergm_out$prior.parameters <- list(prior.mean=prior.mean,prior.cov=priorcov)

  return(BFergm_out)
}

#' @method get_estimates bergm
#' @export
get_estimates.bergm <- function(x, ...){

  K1 <- length(apply(x$Theta,2,median))
  names.bergm.coef <- paste0("theta",1:K1)
  x_MPLE <- ergm::ergmMPLE(formula=x$formula,output="dyadlist")
  design.X <- x_MPLE$predictor[,2+1:K1]
  which.edges <- which(colnames(design.X)=="edges")
  out <- list()
  if(length(which.edges)==0){ #no intercept 'edges'
    out$estimate <- apply(x$Theta,2,median)
    names(out$estimate) <- names.bergm.coef
    mlm1 <- lm(x$Theta ~ 1)
    out$Sigma <- list(sandwich(mlm1) * nrow(x$Theta))
    colnames(out$Sigma[[1]]) <- row.names(out$Sigma[[1]]) <- names(out$estimate)
  }else{
    out$estimate <- apply(x$Theta,2,median)
    names(out$estimate) <- names.bergm.coef
    out$estimate <- out$estimate[-which.edges]
    mlm1 <- lm(x$Theta ~ 1)
    out$Sigma <- list(sandwich(mlm1)[-which.edges,-which.edges] * nrow(x$Theta))
    colnames(out$Sigma[[1]]) <- row.names(out$Sigma[[1]]) <- names(out$estimate)
  }
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "ergm"
  out
}




