#


#' @importFrom stats qnorm dnorm pnorm
#' @describeIn BF S3 method for a named vector 'x'
#' @method BF default
#' @export
BF.default <- function(x,
                       hypothesis = NULL,
                       prior.hyp = NULL,
                       complement = TRUE,
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
                                         prior.hyp = prior.hyp,
                                         complement = complement)

  BF_out$model <- x
  BF_out$call <- match.call()
  BF_out$bayesfactor <- "adjusted fractional Bayes factors using Gaussian approximations"
  BF_out$parameter <- "general parameters"

  BF_out

}

# compute relative meausures (fit or complexity) under a multivariate Gaussian distribution
#' @importFrom mvtnorm dmvnorm pmvnorm
#' @importFrom Matrix rankMatrix
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
    relE <- dmvnorm(rE1,mean=c(meanE),sigma=SigmaE,log=FALSE)
  }
  if(is.null(RrE1) && !is.null(RrO1)){ #only order constraints
    RO1 <- RrO1[,-(K+1)]
    if(!is.matrix(RO1)){
      RO1 <- matrix(RO1,ncol=K)
    }
    qO1 <- nrow(RO1)
    rO1 <- RrO1[,(K+1)]

    if(rankMatrix(RO1)[[1]]==nrow(RO1)){ #RO1 is of full row rank. So use transformation.
      meanO <- c(RO1%*%mean1)
      SigmaO <- RO1%*%Sigma1%*%t(RO1)
      check_vcov(SigmaO)
      relO <- pmvnorm(lower=rO1,upper=Inf,mean=meanO,sigma=SigmaO)[1]
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
    Rr1 <- rbind(RrE1,RrO1)

    if(rankMatrix(Rr1)[[1]] == nrow(Rr1)){

      R1 <- rbind(RE1,RO1)

      #b)
      Tmean1 <- R1 %*% mean1
      TSigma1 <- R1 %*% Sigma1 %*% t(R1)

      # relative meausure for equalities
      relE <- dmvnorm(x=rE1,mean=Tmean1[1:qE1],sigma=matrix(TSigma1[1:qE1,1:qE1],ncol=qE1),log=FALSE)

      # Partitioning equality part and order part
      Tmean1E <- Tmean1[1:qE1]
      Tmean1O <- Tmean1[qE1+1:qO1]

      TSigma1EE <- TSigma1[1:qE1,1:qE1]
      TSigma1OE <- matrix(c(TSigma1[qE1+1:qO1,1:qE1]),nrow=qO1)
      TSigma1OO <- TSigma1[qE1+1:qO1,qE1+1:qO1]

      #conditional location and covariance matrix
      Tmean1OgE <- Tmean1O + TSigma1OE %*% solve(TSigma1EE) %*% (rE1-Tmean1E)
      TSigma1OgE <- TSigma1OO - TSigma1OE %*% solve(TSigma1EE) %*% t(TSigma1OE)

      relO <- pmvnorm(lower=rO1,upper=Inf,mean=c(Tmean1OgE),sigma=TSigma1OgE)[1]

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
#' @importFrom mvtnorm rmvnorm
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
      randomDraws <- rmvnorm(draws2,mean=rep(0,numpara),sigma=diag(numpara))
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

          randomDraws <- rmvnorm(draws2,mean=mean1,sigma=Sigma1)
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




