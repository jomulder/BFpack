#' @importFrom pracma rref
#' @importFrom mvtnorm dmvnorm pmvnorm dmvt pmvt
#' @importFrom Matrix rankMatrix
#' @importFrom MASS ginv
#' @method BF glm
#' @export

source("./R/BFcorrelation.R")
source("./R/BFregression.R") #for make_Rrlist function

BF.glm <- function(x,
                   hypothesis = "exploratory", #instead of using NULL as default, use the actual default
                   prior = NULL,
                   ...){ #I do not see the point of having the dots unless you plan sub-functions of this function
                          #but add them for consistency with BF.lm

  betahat <- x$coefficients
  sigma <- vcov(x)
  K <- length(varnames)
  N <- length(x$residuals)

  varnames <- variable.names(x) #provides the variable names of the object, including intercept
  constraints <- parse_hypothesis(varnames, hypothesis)
  RrList <- make_RrList2(constraints) #make_RrList 1 was not working
  RrE <- do.call(rbind, RrList[[1]])
  RrO <- do.call(rbind, RrList[[2]])

  mean0 <- as.matrix(rep(0,K))
  meanN <- as.matrix(c(betahat))

  prior <- Gaussian_measures(mean0, sigma*N, RrE, RrO)
  posterior <- Gaussian_measures(meanN, sigma, RrE, RrO)


 # parse_hypothesis(constraints) -> list of matrices [function already written]
#  extractor(model) #To be written/copy pasted function to extract means + covariance matrix + N from glm model

 # prior = Gaussian_measures(mean0, sigma1*N, constraints1, constraints2) [Gaussian_measures is the core function in the package everything else is just wrappers]
#  posterior = Gaussian_measures(mean1, sigma1, constraints1, constraints2)

  #From here we just need to compute the output -> can just copy paste from BF.lm (don't forget complement, using function Gaussian_prob_Hc)
#  posterior / prior = BF
#  Gaussian_prob_Hc
#  output-> line 313 (Bayes factors, posterior probabilites)








}


Gaussian_measures <- function(mean1,Sigma1,RrE1,RrO1){
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
      # bain1 <- bain::bain(mean1,Sigma1,RrE1,RrO1,n=10) # choice of n does not matter
      # extract posterior probability (Fit_eq) from bain-object)
      stop("REMINDER. This case should still be implemented.")
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
      TSigma1OE <- TSigma1[qE1+1:qO1,1:qE1]
      TSigma1OO <- TSigma1[qE1+1:qO1,qE1+1:qO1]

      #conditional location and covariance matrix
      Tmean1OgE <- Tmean1O + TSigma1OE %*% solve(TSigma1EE) %*% (rE1-Tmean1E)
      TSigma1OgE <- TSigma1OO - TSigma1OE %*% solve(TSigma1EE) %*% TSigma1OE ##Removed t [transposing] because already correct format?

      relO <- mvtnorm::pmvnorm(lower=rO1,upper=Inf,mean=c(Tmean1OgE),sigma=TSigma1OgE)[1]

    }else{ #use bain for the computation of the probability

      # bain1 <- bain::bain(mean1,Sigma1=Sigma1,RrE1=NULL,RrO1=RO1tilde,n=10) # choice of n does not matter
      # extract posterior probability (Fit_eq) from bain-object)
      stop("REMINDER. This case should still be implemented.")
    }
  }

  return(c(relE,relO))
}
