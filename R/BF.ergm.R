#BF method for ergm class

#' @importFrom sandwich sandwich
#' @importFrom ergm ergmMPLE
#' @importFrom stats as.formula
#' @importFrom Bergm bergm
#' @method BF ergm
#' @export
BF.ergm <- function(x,
                    hypothesis = NULL,
                    prior.hyp.explo = NULL,
                    prior.hyp.conf = NULL,
                    prior.hyp = NULL,
                    complement = TRUE,
                    log = FALSE,
                    cov.prob = .95,
                    ...){

  if(!(cov.prob>0 & cov.prob<1)){
    stop("The argument 'cov.prob' is a coverage probability for the interval estimates that should lie between 0 and 1. The default is 0.95.")
  }

  logIN <- log

  # extract coefficients
  estimate <- coef(x)
  K1 <- length(estimate)
  # get design matrix of pseudo likelihood to construct prior covariance matrix
  nw <- x$network
  form.char <- paste(format(x$formula), collapse = '')
  location_tilde <- regexpr("~",form.char)[1]
  form.new <- as.formula(paste0("nw ~",substr(form.char,start=location_tilde+1,stop=nchar(form.char))))
  x_MPLE <- ergmMPLE(form.new,output="dyadlist")
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
  Bergm.out <- bergm(form.new,prior.mean=rep(0,K1),prior.sigma=priorcov,...)
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
                     prior.hyp.explo = prior.hyp.explo,
                     prior.hyp.conf = prior.hyp.conf,
                     complement = complement,
                     log = logIN)
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
                                         prior.hyp.explo = prior.hyp.explo,
                                         prior.hyp.conf = prior.hyp.conf,
                                         complement = complement,
                                         log = logIN)
  }

  CrI_LB <- (1 - cov.prob)/2
  CrI_UB <- 1 - (1 - cov.prob)/2
  postestimates <- cbind(apply(Bergm.out$Theta,2,mean),
                          apply(Bergm.out$Theta,2,median),
                          apply(Bergm.out$Theta,2,quantile,CrI_LB),
                          apply(Bergm.out$Theta,2,quantile,CrI_UB))
  rownames(postestimates) <- names(estimate)
  colnames(postestimates) <- c("mean","median","2.5%","97.5%")

  BFergm_out$estimates <- postestimates
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

  nw <- x$network
  form.char <- paste(format(x$formula), collapse = '')
  location_tilde <- regexpr("~",form.char)[1]
  form.new <- as.formula(paste0("nw ~",substr(form.char,start=location_tilde+1,stop=nchar(form.char))))

  estimate <- coef(x)
  K1 <- length(estimate)
  x_MPLE <- ergmMPLE(form.new,output="dyadlist")
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
                    prior.hyp.explo = NULL,
                    prior.hyp.conf = NULL,
                    prior.hyp = NULL,
                    complement = TRUE,
                    log = FALSE,
                    ...){

  logIN <- log

  form.char <- paste(format(x$formula), collapse = '')
  location_tilde <- regexpr("~",form.char)[1]
  name.nw <- substr(form.char,start=1,stop=location_tilde-2)
  if(!exists(name.nw)){
    stop(paste0("For an object of class 'bergm', the function 'BF()' only runs if the network data object '",name.nw,
                   "' is also present in the environment."))
  }

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
  x_MPLE <- ergmMPLE(formula=x$formula,output="dyadlist")
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
  Bergm.out <- bergm(x$formula,prior.mean=rep(0,K1),prior.sigma=priorcov,...)
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
                                         prior.hyp.explo = prior.hyp.explo,
                                         prior.hyp.conf = prior.hyp.conf,
                                         complement = complement,
                                         log = logIN)
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
                                         prior.hyp.explo = prior.hyp.explo,
                                         prior.hyp.conf = prior.hyp.conf,
                                         complement = complement,
                                         log = logIN)
  }

  postestimates <- cbind(apply(Bergm.out$Theta,2,mean),
                         apply(Bergm.out$Theta,2,median),
                         apply(Bergm.out$Theta,2,quantile,.025),
                         apply(Bergm.out$Theta,2,quantile,.975))
  rownames(postestimates) <- names(estimate)
  colnames(postestimates) <- c("mean","median","2.5%","97.5%")

  BFergm_out$estimates <- postestimates
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

  form.char <- paste(format(x$formula), collapse = '')
  location_tilde <- regexpr("~",form.char)[1]
  name.nw <- substr(form.char,start=1,stop=location_tilde-2)
  if(!exists(name.nw)){
    stop(paste0("For an object of class 'bergm', the function 'BF()' only runs if the network data object '",name.nw,
                "' is also present in the environment."))
  }

  K1 <- length(apply(x$Theta,2,median))
  names.bergm.coef <- paste0("theta",1:K1)
  x_MPLE <- ergmMPLE(formula=x$formula,output="dyadlist")
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




