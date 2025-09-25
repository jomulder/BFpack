#


#' @importFrom stats qnorm dnorm pnorm
#' @describeIn BF S3 method for a named vector 'x'
#' @method BF default
#' @export
BF.default <- function(x,
                       hypothesis = NULL,
                       prior.hyp.explo = NULL,
                       prior.hyp.conf = NULL,
                       prior.hyp = NULL,
                       complement = TRUE,
                       log = FALSE,
                       cov.prob = .95,
                       Sigma,
                       n,
                       ...){

  #Input is a named mean vector x, covariance matrix and number of observations
  #These are extracted by the relevant method functions from a model object and
  #passed together with the hypothesis and prior to the Gaussian_estimator

  # use Savage-Dickey approximation of the BF
  BF_out <- Savage.Dickey.Gaussian(prior.mean = rep(0, length(x)),
                                        prior.sigma = Sigma * n,
                                        post.mean = x,
                                        post.sigma = Sigma,
                                        hypothesis = hypothesis,
                                        prior.hyp.explo = prior.hyp.explo,
                                        prior.hyp.conf = prior.hyp.conf,
                                        prior.hyp = prior.hyp,
                                        complement = complement,
                                        log = log,
                                   cov.prob=cov.prob)

  BF_out$model <- x
  BF_out$call <- match.call()
  BF_out$bayesfactor <- "adjusted fractional Bayes factors using Gaussian approximations"
  BF_out$parameter <- "general parameters"

  BF_out

}

# extended Savage-Dickey density ratio for multivariate normal prior and posterior
Savage.Dickey.Gaussian <- function(prior.mean,
                                   prior.sigma,
                                   post.mean,
                                   post.sigma,
                                   hypothesis,
                                   prior.hyp.explo,
                                   prior.hyp.conf,
                                   prior.hyp,
                                   complement,
                                   log = FALSE,
                                   cov.prob = .95){

  if(!(cov.prob>0 & cov.prob<1)){
    stop("The argument 'cov.prob' is a coverage probability for the interval estimates that should lie between 0 and 1. The default is 0.95.")
  }
  CrI_LB <- (1 - cov.prob)/2
  CrI_UB <- 1 - (1 - cov.prob)/2

  #prior.mean is a normal prior mean of key parameters
  #prior.sigma is a normal prior covariance matrix of key parameters
  #post.mean is a normal posterior mean of key parameters
  #post.sigma is a normal posterior covariance matrix of key parameters
  #These are passed together with the hypothesis and prior to the Gaussian_estimator
  #the prior of the nuisance parameters under the constrained models is a conditional version of the full prior

  meanN <- post.mean     #the posterior mean is approximated with the estimate
  covmN <- post.sigma    #the posterior covariance matrix is approximated with the error covariance matrix

  names_coef <- names(meanN)
  covm0 <- prior.sigma
  mean0 <- prior.mean # for constrained testing prior mean is relocated to 'boundary of constrained space'

  logIN <- log

  # check proper usage of argument 'prior.hyp.conf' and 'prior.hyp.explo'
  if(!is.null(prior.hyp.conf)){
    prior.hyp <- prior.hyp.conf
  }
  dummy <- 123
  class(dummy) <- "dummy"
  prior.hyp.explo <- process.prior.hyp.explo(prior_hyp_explo = prior.hyp.explo, model=dummy)

  # compute exploratory BFs for each parameter
  relfit <- matrix(c(dnorm(0,mean=meanN,sd=sqrt(diag(covmN)),log=TRUE),
                     pnorm(0,mean=meanN,sd=sqrt(diag(covmN)),log.p=TRUE),
                     pnorm(0,mean=meanN,sd=sqrt(diag(covmN)),log.p=TRUE,lower.tail = FALSE)),ncol=3)

  relcomp <- matrix(c(dnorm(0,mean=mean0,sd=sqrt(diag(covm0)),log=TRUE),
                      pnorm(0,mean=mean0,sd=sqrt(diag(covm0)),log.p=TRUE),
                      pnorm(0,mean=mean0,sd=sqrt(diag(covm0)),log.p=TRUE,lower.tail = FALSE)),ncol=3)
  BFtu_exploratory <- relfit - relcomp
  maxrow <- apply(BFtu_exploratory,1,max)
  norm_BF_explo <- exp(BFtu_exploratory - maxrow %*% t(rep(1,3))) * (rep(1,nrow(relcomp)) %*% t(prior.hyp.explo[[1]]))
  PHP_exploratory <- norm_BF_explo / apply(norm_BF_explo,1,sum)
  # PHP_exploratory <- round(exp(BFtu_exploratory - maxrow %*% t(rep(1,3))) /
  #                            apply(exp(BFtu_exploratory - maxrow %*% t(rep(1,3))),1,sum),3)
  colnames(PHP_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
  row.names(PHP_exploratory) <- row.names(BFtu_exploratory) <- names_coef
  colnames(BFtu_exploratory) <- c("H(=0)","H(<0)","H(>0)")

  # compute posterior estimates
  postestimates <- cbind(meanN,meanN,
                         t(matrix(unlist(lapply(1:length(meanN),function(coef){
                           ub <- qnorm(p=.975)*sqrt(covmN[coef,coef])+meanN[coef]
                           lb <- qnorm(p=.025)*sqrt(covmN[coef,coef])+meanN[coef]
                           return(c(lb,ub))
                         })),nrow=2)),
                         1-pnorm(0,mean=meanN,sd=sqrt(diag(covmN)))
  )
  row.names(postestimates) <- names_coef
  colnames(postestimates) <- c("mean","median",paste0(as.character(round(CrI_LB*100,7)),"%"),
                               paste0(as.character(round(CrI_UB*100,7)),"%"),"Pr(>0)")

  if(logIN == FALSE){
    BFtu_exploratory <- exp(BFtu_exploratory)
  }

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
      #Note that input is a bit strange here, Gaussian_prob_Hc needs fixing
      relcomp <- Gaussian_prob_Hc(mean1 = mean0, Sigma1 = covm0, relmeas = relcomp, RrO = RrO)
      relfit <- Gaussian_prob_Hc(mean1 = meanN, Sigma1 = covmN, relmeas = relfit, RrO = RrO)
    }
    hypothesisshort <- unlist(lapply(1:nrow(relfit),function(h) paste0("H",as.character(h))))
    row.names(relfit) <- row.names(relfit) <- hypothesisshort
    colnames(relcomp) <- c("c_E", "c_0")
    colnames(relfit) <- c("f_E", "f_0")

    # the BF for the complement hypothesis vs Hu needs to be computed.
    BFtu_confirmatory <- c(apply(relfit - relcomp, 1, sum))
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

    maxBFtu <- max(BFtu_confirmatory)
    PHP_confirmatory <- round(exp(BFtu_confirmatory-maxBFtu)*priorprobs /
                                sum(exp(BFtu_confirmatory-maxBFtu)*priorprobs),3)
    BFtable <- cbind(relcomp,relfit,relfit[,1]-relcomp[,1],relfit[,2]-relcomp[,2],
                     apply(relfit,1,sum)-apply(relcomp,1,sum),PHP_confirmatory)
    BFtable[,1:7] <- exp(BFtable[,1:7])
    row.names(BFtable) <- names(PHP_confirmatory)
    colnames(BFtable) <- c("complex=","complex>","fit=","fit>","BF=","BF>","BF","PHP")
    BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)) -
      t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
    diag(BFmatrix_confirmatory) <- log(1)
    # row.names(BFmatrix_confirmatory) <- Hnames
    # colnames(BFmatrix_confirmatory) <- Hnames
    if(nrow(relfit)==length(parse_hyp$original_hypothesis)){
      hypotheses <- parse_hyp$original_hypothesis
    }else{
      hypotheses <- c(parse_hyp$original_hypothesis,"complement")
    }

    if(logIN == FALSE){
      BFtu_confirmatory <- exp(BFtu_confirmatory)
      BFmatrix_confirmatory <- exp(BFmatrix_confirmatory)
    }
  }

  BF_out <- list(
    BFtu_exploratory=BFtu_exploratory,
    PHP_exploratory=PHP_exploratory,
    BFtu_confirmatory=BFtu_confirmatory,
    PHP_confirmatory=PHP_confirmatory,
    BFmatrix_confirmatory=BFmatrix_confirmatory,
    BFtable_confirmatory=BFtable,
    prior.hyp.explo=prior.hyp.explo,
    prior.hyp.conf=priorprobs,
    hypotheses=hypotheses,
    estimates=postestimates,
    model=NULL,
    bayesfactor="Bayes factors using Gaussian approximations",
    parameter="general parameters",
    log = logIN,
    call=NULL)

  class(BF_out) <- "BF"

  BF_out

}

# compute relative meausures (fit or complexity) under a multivariate Gaussian distribution
#' @importFrom mvtnorm dmvnorm pmvnorm
#' @importFrom berryFunctions is.error
Gaussian_measures <- function(mean1,Sigma1,n1=0,RrE1,RrO1,names1=NULL,constraints1=NULL){
  K <- length(mean1)
  relE <- relO <- log(1)
  if(!is.null(RrE1) && is.null(RrO1)){ #only equality constraints
    RE1 <- RrE1[,-(K+1)]
    if(!is.matrix(RE1)){
      RE1 <- matrix(RE1,ncol=K)
    }
    rE1 <- RrE1[,(K+1)]
    qE1 <- nrow(RE1)
    meanE <- RE1%*%mean1
    SigmaE <- RE1%*%Sigma1%*%t(RE1)
    relE <- dmvnorm(rE1,mean=c(meanE),sigma=SigmaE,log=TRUE)
  }
  if(is.null(RrE1) && !is.null(RrO1)){ #only order constraints
    RO1 <- RrO1[,-(K+1)]
    if(!is.matrix(RO1)){
      RO1 <- matrix(RO1,ncol=K)
    }
    qO1 <- nrow(RO1)
    rO1 <- RrO1[,(K+1)]

    if(Rank(RO1)==nrow(RO1)){ #RO1 is of full row rank. So use transformation.
      meanO <- c(RO1%*%mean1)
      SigmaO <- RO1%*%Sigma1%*%t(RO1)
      #check_vcov(SigmaO)
      relO <- log(pmvnorm(lower=rO1,upper=Inf,mean=meanO,sigma=SigmaO)[1])
    }else{ #no linear transformation can be used; pmvt cannot be used. Use bain with a multivariate normal approximation
      names(mean1) <- names1
      mean1vec <- c(mean1)
      names(mean1vec) <- row.names(Sigma1) <- colnames(Sigma1) <- names1
      if(n1>0){ # we need prior measures
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=Sigma1,n=n1)
        relO <- log(bain_res$fit[1,4])
      }else { # we need posterior measures (there is very little information)
        if(!is.error(bain(x=mean1vec,hypothesis=constraints1,Sigma=Sigma1,n=999))){
          bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=Sigma1,n=999) #n not used in computation
          relO <- log(bain_res$fit[1,3])
        }else{
          # simply compute the proportion based on 1e6 draws
          relO <- log(mean(apply(rmvnorm(1e6,mean=mean1vec,sigma=Sigma1) %*% t(RO1) - rep(1,1e6) %*% t(rO1) > 0,1,prod)))
        }
      }

    }
  }
  if(!is.null(RrE1) && !is.null(RrO1)){ #hypothesis with equality and order constraints

    RE1 <- RrE1[,-(K+1)]
    if(!is.matrix(RE1)){
      RE1 <- matrix(RE1,ncol=K)
    }
    rE1 <- RrE1[,(K+1)]
    qE1 <- nrow(RE1)
    RO1 <- RrO1[,-(K+1)]
    if(!is.matrix(RO1)){
      RO1 <- matrix(RO1,ncol=K)
    }
    qO1 <- nrow(RO1)
    rO1 <- RrO1[,(K+1)]
    Rr1 <- rbind(RrE1,RrO1)

    if(Rank(Rr1) == nrow(Rr1)){

      R1 <- rbind(RE1,RO1)

      #b)
      Tmean1 <- R1 %*% mean1
      TSigma1 <- R1 %*% Sigma1 %*% t(R1)

      # relative meausure for equalities
      relE <- dmvnorm(x=rE1,mean=Tmean1[1:qE1],sigma=matrix(TSigma1[1:qE1,1:qE1],ncol=qE1),log=TRUE)

      # Partitioning equality part and order part
      Tmean1E <- Tmean1[1:qE1]
      Tmean1O <- Tmean1[qE1+1:qO1]

      TSigma1EE <- TSigma1[1:qE1,1:qE1]
      TSigma1OE <- matrix(c(TSigma1[qE1+1:qO1,1:qE1]),nrow=qO1)
      TSigma1OO <- TSigma1[qE1+1:qO1,qE1+1:qO1]

      #conditional location and covariance matrix
      Tmean1OgE <- Tmean1O + TSigma1OE %*% solve(TSigma1EE) %*% (rE1-Tmean1E)
      TSigma1OgE <- TSigma1OO - TSigma1OE %*% solve(TSigma1EE) %*% t(TSigma1OE)

      relO <- log(pmvnorm(lower=rO1,upper=Inf,mean=c(Tmean1OgE),sigma=TSigma1OgE)[1])

    }else{ #use bain for the computation of the probability
      names(mean1) <- names1
      if(n1>0){ # we need prior measures
        bain_res <- bain(x=c(mean1),hypothesis=constraints1,Sigma=Sigma1,n=n1)
        relO <- log(bain_res$fit[1,4])
        relE <- log(bain_res$fit[1,2])
      }else { # we need posterior measures (there is very little information)
        bain_res <- bain(x=c(mean1),hypothesis=constraints1,Sigma=Sigma1,n=999) #n not used in computation
        relO <- log(bain_res$fit[1,3])
        relE <- log(bain_res$fit[1,1])
      }
    }
  }

  return(c(relE,relO))
}

# The function computes the probability of an unconstrained draw falling in the complement subspace under a multivariate Gaussian distribution.
#' @importFrom mvtnorm rmvnorm
Gaussian_prob_Hc <- function(mean1,Sigma1,relmeas,RrO){

  numpara <- length(mean1)
  numhyp <- nrow(relmeas)
  which_eq <- relmeas[,1] != log(1)
  if(sum(which_eq)==numhyp){ # Then the complement is equivalent to the unconstrained hypothesis.
    relmeas <- rbind(relmeas,rep(log(1),2))
    rownames(relmeas)[numhyp+1] <- "complement"
  }else{ # So there is at least one hypothesis with only order constraints
    welk <- which(!which_eq)
    if(length(welk)==1){ # There is one hypothesis with only order constraints. Hc is complement of this hypothesis.
      relmeas <- rbind(relmeas,rep(log(1),2))
      relmeas[numhyp+1,2] <- log(1 - exp(relmeas[welk,2]))
      rownames(relmeas)[numhyp+1] <- "complement"
    }else{ # So more than one hypothesis with only order constraints
      # First we check whether ther is an overlap between the order constrained spaces.
      draws2 <- 1e5
      randomDraws <- rmvnorm(draws2,mean=rep(0,numpara),sigma=diag(numpara))
      #get draws that satisfy the constraints of the separate order constrained hypotheses
      checksOC <- lapply(welk,function(h){
        Rorder <- matrix(RrO[[h]][,-(1+numpara)],nrow=nrow(RrO[[h]]))
        # if(ncol(Rorder)==1){
        #   Rorder <- t(Rorder)
        # }
        rorder <- as.matrix(RrO[[h]][,1+numpara])
        apply(randomDraws%*%t(Rorder) > rep(1,draws2)%*%t(rorder),1,prod)
      })
      checkOCplus <- Reduce("+",checksOC)

      if(sum(checkOCplus > 0) < draws2){ #then the joint order constrained hypotheses do not completely cover the parameter space.
        if(sum(checkOCplus>1)==0){ # then order constrained spaces are nonoverlapping
          relmeas <- rbind(relmeas,rep(log(1),2))
          relmeas[numhyp+1,2] <- log(1 - sum(exp(relmeas[welk,2])))
          rownames(relmeas)[numhyp+1] <- "complement"
        }else{ #the order constrained subspaces at least partly overlap

          # funtion below gives a rough estimate of the posterior probability under Hc
          # a bain type of algorithm would be better of course. but for now this is ok.

          randomDraws <- rmvnorm(draws2,mean=mean1,sigma=Sigma1)
          checksOCpost <- lapply(welk,function(h){
            Rorder <- matrix(RrO[[h]][,-(1+numpara)],nrow=nrow(RrO[[h]]))
            # if(ncol(Rorder)==1){
            #   Rorder <- t(Rorder)
            # }
            rorder <- as.matrix(RrO[[h]][,1+numpara])
            apply(randomDraws%*%t(Rorder) > rep(1,draws2)%*%t(rorder),1,prod)
          })
          relmeas <- rbind(relmeas,rep(log(1),2))
          relmeas[numhyp+1,2] <- log(sum(Reduce("+",checksOCpost) == 0) / draws2)
          rownames(relmeas)[numhyp+1] <- "complement"
        }
      }
    }
  }

  return(relmeas)
}




