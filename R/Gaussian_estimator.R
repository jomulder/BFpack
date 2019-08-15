##Internal estimation function for the methods for classes glm, lavaan, coxph, rem, rem.dyad and glmerMod

#' @importFrom stats qnorm
BF_Gaussian <- function(meanN,
                       covmN,
                       n,
                       hypothesis,
                       prior){

  #Input is a named mean vector, covariance matrix and number of observations
  #These are extracted by the relevant method functions from a model object and
  #passed together with the hypothesis and prior to the Gaussian_estimator

  names_coef <- names(meanN)
  covm0 <- covmN * n
  mean0 <- as.matrix(rep(0, length(names_coef)))

  # compute exploratory BFs for each parameter
  relfit <- matrix(c(dnorm(0,mean=meanN,sd=sqrt(diag(covmN))), #[Anton] Are these relfit/relcomp computations general or specific to correlations? [Joris] This is general. So it also works for these parameters.
                     pnorm(0,mean=meanN,sd=sqrt(diag(covmN))),
                     1-pnorm(0,mean=meanN,sd=sqrt(diag(covmN)))),ncol=3)

  relcomp <- matrix(c(dnorm(0,mean=mean0,sd=sqrt(diag(covm0))),
                      pnorm(0,mean=mean0,sd=sqrt(diag(covm0))),
                      1-pnorm(0,mean=mean0,sd=sqrt(diag(covm0)))),ncol=3)
  BFtu_exploratory <- relfit / relcomp
  PHP_exploratory <- round(BFtu_exploratory / apply(BFtu_exploratory,1,sum),3)
  colnames(PHP_exploratory) <- c("p(=0)","Pr(<0)","Pr(>0)")
  row.names(PHP_exploratory) <- names_coef

  # compute posterior estimates
  postestimates <- cbind(meanN,meanN,
                         t(matrix(unlist(lapply(1:length(meanN),function(coef){
                           ub <- qnorm(p=.975)*sqrt(covmN[coef,coef])+meanN[coef]
                           lb <- qnorm(p=.025)*sqrt(covmN[coef,coef])+meanN[coef]
                           return(c(lb,ub))
                         })),nrow=2))
  )
  row.names(postestimates) <- names_coef
  colnames(postestimates) <- c("mean","median","2.5%","97.5%")

  if(is.null(hypothesis)){
    BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relfit <-
      relcomp <- BFtable <- hypotheses <- priorprobs <- NULL
  }else{
    # confirmatory tests based on input constraints
    parse_hyp <- parse_hypothesis(names_coef,hypothesis)

    #create coefficient with equality and order constraints
    RrList <- make_RrList2(parse_hyp)
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]

    # RrStack is used to check conflicting constraints, and for the default prior location
    RrStack <- rbind(do.call(rbind,RrE),do.call(rbind,RrO))
    RrStack <- interval_RrStack(RrStack)
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
      nonzero <- RrStack[,K+1]!=0
      if(max(nonzero)>0){
        row1 <- max(which(nonzero))
        if(sum(abs(rref_ei[row1,1:K]))==0){
          stop("No common boundary point for prior location. Conflicting constraints.")
        }
      }
    }
    # default prior location
    mean0 <- ginv(RStack)%*%rStack

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

    # get relative fit and complexity of complement hypothesis
    relcomp <- Gaussian_prob_Hc(mean1 = mean0, Sigma1 = covm0, relmeas = relcomp, RrO = RrO) #Note that input is a bit strange here, Gaussian_prob_Hc needs fixing
    relfit <- Gaussian_prob_Hc(mean1 = meanN, Sigma1 = covmN, relmeas = relfit, RrO = RrO)

    # Hnames <- c(unlist(lapply(1:numhyp,function(h){paste0("H",as.character(h))})),"Hc")
    # row.names(relcomp) <- row.names(relfit) <- Hnames
    Hnames <- row.names(relcomp)
    colnames(relcomp) <- c("c_E", "c_0")
    colnames(relfit) <- c("f_E", "f_0")

    # the BF for the complement hypothesis vs Hu needs to be computed.
    BFtu_confirmatory <- c(apply(relfit / relcomp, 1, prod))
    # Check input of prior probabilies
    if(is.null(prior)){
      priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
    }else{
      if(!is.numeric(prior) || length(prior)!=length(BFtu_confirmatory)){
        warning(paste0("Argument 'prior' should be numeric and of length ",as.character(length(BFtu_confirmatory)),". Equal prior probabilities are used."))
        priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
      }else{
        priorprobs <- prior
      }
    }

    PHP_confirmatory <- round(BFtu_confirmatory*priorprobs / sum(BFtu_confirmatory*priorprobs),3)
    BFtable <- cbind(relcomp,relfit,relfit[,1]/relcomp[,1],relfit[,2]/relcomp[,2],
                     apply(relfit,1,prod)/apply(relcomp,1,prod),PHP_confirmatory)
    row.names(BFtable) <- names(PHP_confirmatory)
    colnames(BFtable) <- c("comp_E","comp_O","fit_E","fit_O","BF_E","BF_O","BF","PHP")
    BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory))/
      t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
    row.names(BFmatrix_confirmatory) <- Hnames
    colnames(BFmatrix_confirmatory) <- Hnames
    hypotheses <- Hnames
  }

    BF_out <- list(
      BFtu_exploratory=BFtu_exploratory,
      PHP_exploratory=PHP_exploratory,
      BFtu_confirmatory=BFtu_confirmatory,
      PHP_confirmatory=PHP_confirmatory,
      BFmatrix_confirmatory=BFmatrix_confirmatory,
      BFtable_confirmatory=BFtable,
      prior=priorprobs,
      hypotheses=hypotheses,
      estimates=postestimates)

    class(BF_out) <- "BF"

    BF_out

}



# compute relative meausures (fit or complexity) under a multivariate Gaussian distribution
Gaussian_measures <- function(mean1,Sigma1,n1=0,RrE1,RrO1,names1=NULL,constraints1=NULL){
  K <- length(mean1)
  relE <- relO <- 1
  if(!is.null(RrE1) && is.null(RrO1)){ #only equality constraints
    RE1 <- RrE1[,-(K+1)]
    if(!is.matrix(RE1)){
      RE1 <- matrix(RE1,ncol=K)
    }
    rE1 <- RrE1[,(K+1)]
    qE1 <- nrow(RE1)
    meanE <- RE1%*%mean1
    SigmaE <- RE1%*%Sigma1%*%t(RE1)
    relE <- mvtnorm::dmvnorm(rE1,mean=c(meanE),sigma=SigmaE,log=FALSE)
  }
  if(is.null(RrE1) && !is.null(RrO1)){ #only order constraints
    RO1 <- RrO1[,-(K+1)]
    if(!is.matrix(RO1)){
      RO1 <- matrix(RO1,ncol=K)
    }
    qO1 <- nrow(RO1)
    rO1 <- RrO1[,(K+1)]

    if(Matrix::rankMatrix(RO1)[[1]]==nrow(RO1)){ #RO1 is of full row rank. So use transformation.
      meanO <- c(RO1%*%mean1)
      SigmaO <- RO1%*%Sigma1%*%t(RO1)
      relO <- mvtnorm::pmvnorm(lower=rO1,upper=Inf,mean=meanO,sigma=SigmaO)[1]
    }else{ #no linear transformation can be used; pmvt cannot be used. Use bain with a multivariate normal approximation
      names(mean1) <- names1
      if(n1>0){ # we need prior measures
        mean1vec <- c(mean1)
        names(mean1vec) <- names1
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=Sigma1,n=n1)
        relO <- bain_res$fit[1,4]
      }else { # we need posterior measures (there is very little information)
        mean1vec <- c(mean1)
        names(mean1vec) <- names1
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=Sigma1,n=999) #n not used in computation
        relO <- bain_res$fit[1,3]
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

    if(Matrix::rankMatrix(RrO1)[[1]] == nrow(RrO1)){

      R1 <- rbind(RE1,RO1)

      #b)
      Tmean1 <- R1 %*% mean1
      TSigma1 <- R1 %*% Sigma1 %*% t(R1)

      # relative meausure for equalities
      relE <- mvtnorm::dmvnorm(x=rE1,mean=Tmean1[1:qE1],sigma=matrix(TSigma1[1:qE1,1:qE1],ncol=qE1),log=FALSE)

      # Partitioning equality part and order part
      Tmean1E <- Tmean1[1:qE1]
      Tmean1O <- Tmean1[qE1+1:qO1]

      TSigma1EE <- TSigma1[1:qE1,1:qE1]
      TSigma1OE <- matrix(c(TSigma1[qE1+1:qO1,1:qE1]),nrow=qO1)
      TSigma1OO <- TSigma1[qE1+1:qO1,qE1+1:qO1]

      #conditional location and covariance matrix
      Tmean1OgE <- Tmean1O + TSigma1OE %*% solve(TSigma1EE) %*% (rE1-Tmean1E)
      TSigma1OgE <- TSigma1OO - TSigma1OE %*% solve(TSigma1EE) %*% t(TSigma1OE)

      relO <- mvtnorm::pmvnorm(lower=rO1,upper=Inf,mean=c(Tmean1OgE),sigma=TSigma1OgE)[1]

    }else{ #use bain for the computation of the probability
      names(mean1) <- names1
      if(n1>0){ # we need prior measures
        bain_res <- bain(x=c(mean1),hypothesis=constraints1,Sigma=Sigma1,n=n1)
        relO <- bain_res$fit[1,4]
        relE <- bain_res$fit[1,2]
      }else { # we need posterior measures (there is very little information)
        bain_res <- bain(x=c(mean1),hypothesis=constraints1,Sigma=Sigma1,n=999) #n not used in computation
        relO <- bain_res$fit[1,3]
        relE <- bain_res$fit[1,1]
      }
    }
  }

  return(c(relE,relO))
}

# The function computes the probability of an unconstrained draw falling in the complement subspace under a multivariate Gaussian distribution.
Gaussian_prob_Hc <- function(mean1,Sigma1,relmeas,RrO){

  numpara <- length(mean1)
  numhyp <- nrow(relmeas)
  which_eq <- relmeas[,1] != 1
  if(sum(which_eq)==numhyp){ # Then the complement is equivalent to the unconstrained hypothesis.
    relmeas <- rbind(relmeas,rep(1,2))
    rownames(relmeas)[numhyp+1] <- "complement"
  }else{ # So there is at least one hypothesis with only order constraints
    welk <- which(!which_eq)
    if(length(welk)==1){ # There is one hypothesis with only order constraints. Hc is complement of this hypothesis.
      relmeas <- rbind(relmeas,rep(1,2))
      relmeas[numhyp+1,2] <- 1 - relmeas[welk,2]
      rownames(relmeas)[numhyp+1] <- "complement"
    }else{ # So more than one hypothesis with only order constraints
      # First we check whether ther is an overlap between the order constrained spaces.
      draws2 <- 1e4
      randomDraws <- mvtnorm::rmvnorm(draws2,mean=rep(0,numpara),sigma=diag(numpara))
      #get draws that satisfy the constraints of the separate order constrained hypotheses
      checksOC <- lapply(welk,function(h){
        Rorder <- as.matrix(RrO[[h]][,-(1+numpara)])
        if(ncol(Rorder)==1){
          Rorder <- t(Rorder)
        }
        rorder <- as.matrix(RrO[[h]][,1+numpara])
        apply(randomDraws%*%t(Rorder) > rep(1,draws2)%*%t(rorder),1,prod)
      })
      checkOCplus <- Reduce("+",checksOC)

      if(sum(checkOCplus > 0) < draws2){ #then the joint order constrained hypotheses do not completely cover the parameter space.
        if(sum(checkOCplus>1)==0){ # then order constrained spaces are nonoverlapping
          relmeas <- rbind(relmeas,rep(1,2))
          relmeas[numhyp+1,2] <- 1 - sum(relmeas[welk,2])
          rownames(relmeas)[numhyp+1] <- "complement"
        }else{ #the order constrained subspaces at least partly overlap

          # funtion below gives a rough estimate of the posterior probability under Hc
          # a bain type of algorithm would be better of course. but for now this is ok.

          randomDraws <- mvtnorm::rmvnorm(draws2,mean=mean1,sigma=Sigma1)
          checksOCpost <- lapply(welk,function(h){
            Rorder <- as.matrix(RrO[[h]][,-(1+numpara)])
            if(ncol(Rorder)==1){
              Rorder <- t(Rorder)
            }
            rorder <- as.matrix(RrO[[h]][,1+numpara])
            apply(randomDraws%*%t(Rorder) > rep(1,draws2)%*%t(rorder),1,prod)
          })
          relmeas <- rbind(relmeas,rep(1,2))
          relmeas[numhyp+1,2] <- sum(Reduce("+",checksOCpost) == 0) / draws2
          rownames(relmeas)[numhyp+1] <- "complement"
        }
      }
    }
  }

  return(relmeas)
}




