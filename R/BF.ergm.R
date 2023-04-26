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

  # first check if effect names in hypothesis argument correspond with names in x
  coef_names <- names(get_estimates(x)$estimate)
  if(!is.null(hypothesis)){
    test0 <- parse_hypothesis(coef_names,hypothesis)
  }

  # extract coefficients
  estimate <- coef(x)
  K1 <- length(estimate)
  # get design matrix of pseudo likelihood
  x_MPLE <- ergm::ergmMPLE(formula=x$formula,output="dyadlist")
  Xdelta <- as.matrix(x_MPLE$predictor[,2+1:K1])
  priorcov <- solve(t(Xdelta)%*%Xdelta) * nrow(Xdelta)
  Bergm.out <- Bergm::bergm(x$formula,prior.mean=rep(0,K1),prior.sigma=priorcov,...)
  #get robust estimates for the Gaussian mean and covariance matrix
  post.mean <- apply(Bergm.out$Theta,2,median)
  names(post.mean) <- names(estimate)
  #get robust estimate of posterior covariance matrix
  mlm1 <- lm(Bergm.out$Theta ~ 1)
  post.Sigma <- sandwich(mlm1) * nrow(Bergm.out$Theta)

  # use Savage-Dickey approximation of the BF
  BFergm_out <- Savage.Dickey.Gaussian(prior.mean = rep(0,K1),
                     prior.sigma = priorcov,
                     post.mean = post.mean,
                     post.sigma = post.Sigma,
                     hypothesis = hypothesis,
                     prior.hyp = prior.hyp,
                     complement = complement)
  approx_estimates <- BFergm_out$estimates
  #update table of estimates with Bayesian estimates
  BFergm_out$estimates <- cbind(apply(Bergm.out$Theta,2,mean),
                                apply(Bergm.out$Theta,2,median),
                                apply(Bergm.out$Theta,2,sd),
                                apply(Bergm.out$Theta,2,quantile,.025),
                                apply(Bergm.out$Theta,2,quantile,.975))
  colnames(BFergm_out$estimates) <- colnames(approx_estimates)
  row.names(BFergm_out$estimates) <- colnames(Xdelta)
  BFergm_out$model <- x
  BFergm_out$call <- match.call()
  BFergm_out$bayesfactor <- "Bayes factors based on unit-information priors and Gaussian approximations"
  BFergm_out$parameter <- "ERGM coefficients"
  BFergm_out$model_update <- Bergm.out

  return(BFergm_out)
}

#' @method get_estimates ergm
#' @export
get_estimates.ergm <- function(x, ...){
  out <- list()
  out$estimate <- coef(x)
  out$Sigma <- list(vcov(x))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "ergm"
  out
}


#' @importFrom ergm ergmMPLE
#' @method BF bergm
#' @export
BF.bergm <- function(x,
                    hypothesis = NULL,
                    prior.hyp = NULL,
                    complement = TRUE,
                    ...){

  # first check if effect names in hypothesis argument correspond with names in x
  coef_names <- names(get_estimates(x)$estimate)
  if(!is.null(hypothesis)){
    test0 <- parse_hypothesis(coef_names,hypothesis)
  }

  K1 <- ncol(x$Theta)
  # get design matrix of pseudo likelihood
  x_MPLE <- ergm::ergmMPLE(formula=x$formula,output="dyadlist")
  Xdelta <- as.matrix(x_MPLE$predictor[,2+1:K1])
  priorcov <- solve(t(Xdelta)%*%Xdelta) * nrow(Xdelta)
  Bergm.out <- Bergm::bergm(x$formula,prior.mean=rep(0,K1),prior.sigma=priorcov,...)
  #get robust estimates for the Gaussian mean and covariance matrix
  post.mean <- apply(Bergm.out$Theta,2,median)
  names(post.mean) <- coef_names
  #get robust estimate of posterior covariance matrix
  mlm1 <- lm(Bergm.out$Theta ~ 1)
  post.Sigma <- sandwich(mlm1) * nrow(Bergm.out$Theta)

  # use Savage-Dickey approximation of the BF
  BFergm_out <- Savage.Dickey.Gaussian(prior.mean = rep(0,K1),
                                     prior.sigma = priorcov,
                                     post.mean = post.mean,
                                     post.sigma = post.Sigma,
                                     hypothesis = hypothesis,
                                     prior.hyp = prior.hyp,
                                     complement = complement)
  approx_estimates <- BFergm_out$estimates
  #update table of estimates with Bayesian estimates
  BFergm_out$estimates <- cbind(apply(Bergm.out$Theta,2,mean),
                                apply(Bergm.out$Theta,2,median),
                                apply(Bergm.out$Theta,2,sd),
                                apply(Bergm.out$Theta,2,quantile,.025),
                                apply(Bergm.out$Theta,2,quantile,.975))
  colnames(BFergm_out$estimates) <- colnames(approx_estimates)
  row.names(BFergm_out$estimates) <- coef_names
  BFergm_out$model <- x
  BFergm_out$call <- match.call()
  BFergm_out$bayesfactor <- "Bayes factors based on unit-information priors and Gaussian approximations"
  BFergm_out$parameter <- "ERGM coefficients"
  BFergm_out$model_update <- Bergm.out

  return(BFergm_out)
}

#' @method get_estimates bergm
#' @export
get_estimates.bergm <- function(x, ...){
  out <- list()
  out$estimate <- apply(x$Theta,2,median)
  K1 <- length(out$estimate)
  names(out$estimate) <- paste0("theta",1:K1)
  mlm1 <- lm(x$Theta ~ 1)
  out$Sigma <- list(sandwich(mlm1) * nrow(x$Theta))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "ergm"
  out
}

# extended Savage-Dickey density ratio for multivariate normal prior and posterior
Savage.Dickey.Gaussian <- function(prior.mean,
                                     prior.sigma,
                                     post.mean,
                                     post.sigma,
                                     hypothesis,
                                     prior.hyp,
                                     complement){

  #prior.mean is a normal prior mean of key parameters
  #prior.sigma is a normal prior covariance matrix of key parameters
  #post.mean is a normal posterior mean of key parameters
  #post.sigma is a normal posterior covariance matrix of key parameters
  #These are passed together with the hypothesis and prior to the Gaussian_estimator

  meanN <- post.mean     #the posterior mean is approximated with the estimate
  covmN <- post.sigma    #the posterior covariance matrix is approximated with the error covariance matrix

  names_coef <- names(meanN)
  covm0 <- prior.sigma
  mean0 <- prior.mean # for constrained testing prior mean is relocated to 'boundary of constrained space'

  # compute exploratory BFs for each parameter
  relfit <- matrix(c(dnorm(0,mean=meanN,sd=sqrt(diag(covmN))),
                     pnorm(0,mean=meanN,sd=sqrt(diag(covmN))),
                     1-pnorm(0,mean=meanN,sd=sqrt(diag(covmN)))),ncol=3)

  relcomp <- matrix(c(dnorm(0,mean=mean0,sd=sqrt(diag(covm0))),
                      pnorm(0,mean=mean0,sd=sqrt(diag(covm0))),
                      1-pnorm(0,mean=mean0,sd=sqrt(diag(covm0)))),ncol=3)
  BFtu_exploratory <- relfit / relcomp
  PHP_exploratory <- round(BFtu_exploratory / apply(BFtu_exploratory,1,sum),3)
  colnames(PHP_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
  row.names(PHP_exploratory) <- names_coef

  # compute posterior estimates
  postestimates <- cbind(meanN,meanN,sqrt(diag(covmN)),
                         t(matrix(unlist(lapply(1:length(meanN),function(coef){
                           ub <- qnorm(p=.975)*sqrt(covmN[coef,coef])+meanN[coef]
                           lb <- qnorm(p=.025)*sqrt(covmN[coef,coef])+meanN[coef]
                           return(c(lb,ub))
                         })),nrow=2))
  )
  row.names(postestimates) <- names_coef
  colnames(postestimates) <- c("mean","median","sd","2.5%","97.5%")

  if(is.null(hypothesis)){
    BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relfit <-
      relcomp <- BFtable <- hypotheses <- priorprobs <- NULL
  }else{
    # confirmatory tests based on input constraints
    parse_hyp <- parse_hypothesis(names_coef,hypothesis)
    parse_hyp$hyp_mat <- do.call(rbind, parse_hyp$hyp_mat)
    #create coefficient with equality and order constraints
    RrList <- make_RrList2(parse_hyp)
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]

    # RrStack is used to check conflicting constraints, and for the default prior location
    if(length(RrE)==1){
      RrStack <- rbind(do.call(rbind,RrE),do.call(rbind,RrO))
      RrStack <- interval_RrStack(RrStack)
    }else{
      RrStack_list <- lapply(1:length(RrE),function(h){
        interval_RrStack(rbind(RrE[[h]],RrO[[h]]))
      })
      RrStack <- do.call(rbind,RrStack_list)
    }
    K <- length(meanN)
    if(nrow(RrStack)>1){
      RStack <- RrStack[,-(K+1)]
      rStack <- RrStack[,(K+1)]
    }else{
      RStack <- matrix(RrStack[,-(K+1)],nrow=1)
      rStack <- RrStack[,(K+1)]
    }

    # check if a common boundary exists for prior location under all constrained hypotheses
    if(nrow(RrStack) > 1){
      rref_ei <- rref(RrStack)
      nonzero <- rref_ei[,K+1]!=0
      if(max(nonzero)>0){
        row1 <- max(which(nonzero))
        if(sum(abs(rref_ei[row1,1:K]))==0){
          stop("No common boundary point for prior location. Conflicting constraints.")
        }
      }
      #determine fraction via number of independent rows (constraints)
      if(is.matrix(rref_ei[,-(K+1)])){
        numindep <- sum(apply(abs(rref_ei[,-(K+1)]),1,sum)!=0)
      }else{
        numindep <- sum(apply(abs(as.matrix(rref_ei[,-(K+1)])),1,sum)!=0)
      }
    } else {
      numindep <- 1
    }

    #get relative fit and complexity of hypotheses
    numhyp <- length(RrE)

    relcomp <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Gaussian_measures(mean1 = mean0, Sigma1 = covm0, RrE1 = RrE[[h]], RrO1 = RrO[[h]],
                        names1=names_coef,constraints1=parse_hyp$original_hypothesis[h])
    })),nrow=2))
    relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Gaussian_measures(mean1 = meanN, Sigma1 = covmN, RrE1 = RrE[[h]], RrO1 = RrO[[h]],
                        names1=names_coef,constraints1=parse_hyp$original_hypothesis[h])
    })),nrow=2))
    row.names(relfit) <- row.names(relcomp) <- parse_hyp$original_hypothesis

    if(complement == TRUE){
      # get relative fit and complexity of complement hypothesis
      relcomp <- Gaussian_prob_Hc(mean1 = mean0, Sigma1 = covm0, relmeas = relcomp, RrO = RrO) #Note that input is a bit strange here, Gaussian_prob_Hc needs fixing
      relfit <- Gaussian_prob_Hc(mean1 = meanN, Sigma1 = covmN, relmeas = relfit, RrO = RrO)
    }
    hypothesisshort <- unlist(lapply(1:nrow(relfit),function(h) paste0("H",as.character(h))))
    row.names(relfit) <- row.names(relfit) <- hypothesisshort
    colnames(relcomp) <- c("c_E", "c_0")
    colnames(relfit) <- c("f_E", "f_0")

    # the BF for the complement hypothesis vs Hu needs to be computed.
    BFtu_confirmatory <- c(apply(relfit / relcomp, 1, prod))
    # Check input of prior probabilies
    if(is.null(prior.hyp)){
      priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
    }else{
      if(!is.numeric(prior.hyp) || length(prior.hyp)!=length(BFtu_confirmatory)){
        warning(paste0("Argument 'prior.hyp' should be numeric and of length ",as.character(length(BFtu_confirmatory)),". Equal prior probabilities are used."))
        priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
      }else{
        priorprobs <- prior.hyp
      }
    }
    names(priorprobs) <- hypothesisshort

    PHP_confirmatory <- round(BFtu_confirmatory*priorprobs / sum(BFtu_confirmatory*priorprobs),3)
    BFtable <- cbind(relcomp,relfit,relfit[,1]/relcomp[,1],relfit[,2]/relcomp[,2],
                     apply(relfit,1,prod)/apply(relcomp,1,prod),PHP_confirmatory)
    row.names(BFtable) <- names(PHP_confirmatory)
    colnames(BFtable) <- c("complex=","complex>","fit=","fit>","BF=","BF>","BF","PHP")
    BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory))/
      t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
    diag(BFmatrix_confirmatory) <- 1
    # row.names(BFmatrix_confirmatory) <- Hnames
    # colnames(BFmatrix_confirmatory) <- Hnames
    if(nrow(relfit)==length(parse_hyp$original_hypothesis)){
      hypotheses <- parse_hyp$original_hypothesis
    }else{
      hypotheses <- c(parse_hyp$original_hypothesis,"complement")
    }
  }

  BF_out <- list(
    BFtu_exploratory=BFtu_exploratory,
    PHP_exploratory=PHP_exploratory,
    BFtu_confirmatory=BFtu_confirmatory,
    PHP_confirmatory=PHP_confirmatory,
    BFmatrix_confirmatory=BFmatrix_confirmatory,
    BFtable_confirmatory=BFtable,
    prior.hyp=priorprobs,
    hypotheses=hypotheses,
    estimates=postestimates,
    model=NULL,
    bayesfactor="Bayes factors using Gaussian approximations",
    parameter="general parameters",
    call=NULL)

  class(BF_out) <- "BF"

  BF_out

}




