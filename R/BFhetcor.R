


#' @importFrom MCMCpack rinvgamma
#' @method BF hetcor
#' @export
BF.hetcor <- function(x,
                       hypothesis = NULL,
                       prior = NULL,
                       ...){
  if(is.null(hypothesis)){
    constraints <- "exploratory"
  } else {
    constraints <- hypothesis
  }
  if(is.null(prior)){
    priorprob <- "default"
  } else {
    priorprob <- prior
  }

  P <- nrow(x$std.errors)
  numcorr <- P*(P-1)/2
  estimates <- x$correlations[lower.tri(diag(P))]
  stderr <- x$std.errors[lower.tri(diag(P))]
  errcov <- diag(stderr**2)

  corr_names <- names(get_estimates(x$correlations)$estimate)
  matrix_names <- matrix(corr_names,nrow=P)
  # equal correlations are at the opposite side of the vector
  corr_names_lower <- matrix_names[lower.tri(matrix_names)]
  corr_names_all <- c(matrix_names[lower.tri(matrix_names)],
                  t(matrix_names)[lower.tri(matrix_names)])

  #exploratory BF testing
  relfit <- matrix(c(dnorm(0,mean=estimates,sd=sqrt(diag(errcov))),
                     pnorm(0,mean=estimates,sd=sqrt(diag(errcov))),
                     1-pnorm(0,mean=estimates,sd=sqrt(diag(errcov)))),ncol=3)
  relcomp <- matrix(c(dbeta(rep(.5,numcorr),shape1=P/2,shape2=P/2)*.5,
                      rep(.5,2*numcorr)),ncol=3)
  row.names(relfit) <- row.names(relcomp) <- corr_names_lower

  BFtu_exploratory <- relfit / relcomp
  colnames(BFtu_exploratory) <- colnames(BFtu_exploratory) <-  c("Pr(=0)","Pr(<0)","Pr(>0)")
  PHP_exploratory <- BFtu_exploratory / apply(BFtu_exploratory,1,sum)

  #confirmatory BF testing
  if(constraints!="exploratory"){

    parse_hyp <- parse_hypothesis(corr_names_all,constraints)
    #combine equivalent correlations, e.g., cor(Y1,Y2)=corr(Y2,Y1).
    parse_hyp$hyp_mat <- cbind(parse_hyp$hyp_mat[,1:numcorr] + parse_hyp$hyp_mat[,numcorr+1:numcorr],
            parse_hyp$hyp_mat[,numcorr*2+1])
    #create coefficient with equality and order constraints
    RrList <- make_RrList2(parse_hyp)
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]

    numhyp <- length(RrE)
    relcomp <- t(matrix(unlist(lapply(1:numhyp,function(h){
      jointuniform_measures(P,numcorr,1,RrE[[h]],RrO[[h]],Fisher=0)
    })),nrow=2))
    relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Gaussian_measures(estimates,errcov,RrE1=RrE[[h]],RrO1=RrO[[h]],names1=corr_names_lower,
                        constraints1=parse_hyp$original_hypothesis[h])
    })),nrow=2))
    row.names(relcomp) <- parse_hyp$original_hypothesis
    row.names(relfit) <- parse_hyp$original_hypothesis
    # evaluation of complement hypothesis
    relfit <- Gaussian_prob_Hc(estimates,errcov,relfit,RrO)
    relcomp <- jointuniform_prob_Hc(P,numcorr,1,relcomp,RrO)

    colnames(relcomp) <- c("c_E","c_O")
    colnames(relfit) <- c("f_E","f_O")
    # computation of exploratory BFs and PHPs
    # the BF for the complement hypothesis vs Hu needs to be computed.
    BFtu_confirmatory <- c(apply(relfit / relcomp, 1, prod))
    # Check input of prior probabilies
    if(!(priorprob == "default" || (length(priorprob)==nrow(relfit) && min(priorprob)>0) )){
      stop("'probprob' must be a vector of positive values or set to 'default'.")
    }
    # Change prior probs in case of default setting
    if(priorprob=="default"){
      priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
    }else{
      priorprobs <- priorprobs/sum(priorprobs)
    }
    names(priorprobs) <- names(BFtu_confirmatory)
    PHP_confirmatory <- BFtu_confirmatory*priorprobs / sum(BFtu_confirmatory*priorprobs)
    # names(PHP_confirmatory) <- unlist(lapply(1:length(parse_hyp$original_hypothesis),function(hyp){
    #   paste0("Pr(",parse_hyp$original_hypothesis,")")
    # }))
    BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory))/
      t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
    row.names(BFmatrix_confirmatory) <- colnames(BFmatrix_confirmatory) <- names(BFtu_confirmatory)
  }else{
    BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relfit <-
      relcomp <- NULL
  }

  # Store in output
  BFlm_out <- list(
    BFtu_exploratory=BFtu_exploratory,
    PHP_exploratory=PHP_exploratory,
    BFtu_confirmatory=BFtu_confirmatory,
    PHP_confirmatory=PHP_confirmatory,
    BFmatrix_confirmatory=BFmatrix_confirmatory,
    relative_fit=relfit,
    relative_complexity=relcomp,
    model=x,
    hypothesis=hypothesis,
    prior=prior)

  class(BFlm_out) <- "BF"

  return(BFlm_out)

}




# compute relative meausures (fit or complexity) under a multivariate Gaussian distribution
Gaussian_measures <- function(mean1,Sigma1,n1=0,RrE1,RrO1,names1=NULL,constraints1=NULL){
  K <- length(mean1)
  relE <- relO <- 1
  if(!is.null(RrE1) && is.null(RrO1)){ #only equality constraints
    RE1 <- RrE1[,-(K+1)]
    if(!is.matrix(RE1)){
      RE1 <- t(as.matrix(RE1))
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
      RO1 <- t(as.matrix(RO1))
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
        bain_res <- bain(x=c(mean1),hypothesis=constraints1,Sigma=Sigma1,n=n1)
        relO <- bain_res$fit[1,4]
      }else { # we need posterior measures (there is very little information)
        bain_res <- bain(x=c(mean1),hypothesis=constraints1,Sigma=Sigma1,n=999) #n not used in computation
        relO <- bain_res$fit[1,3]
      }
    }
  }
  if(!is.null(RrE1) && !is.null(RrO1)){ #hypothesis with equality and order constraints

    RE1 <- RrE1[,-(K+1)]
    if(!is.matrix(RE1)){
      RE1 <- t(as.matrix(RE1))
    }
    rE1 <- RrE1[,(K+1)]
    qE1 <- nrow(RE1)
    RO1 <- RrO1[,-(K+1)]
    if(!is.matrix(RO1)){
      RO1 <- t(as.matrix(RO1))
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
  #relmeas <- relmeas[1:numhyp,]
  which_eq <- relmeas[,1] != 1
  if(sum(which_eq)==numhyp){ # Then the complement is equivalent to the unconstrained hypothesis.
    relmeas <- rbind(relmeas,rep(1,2))
    rownames(relmeas)[relmeas+1] <- "complement"
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



