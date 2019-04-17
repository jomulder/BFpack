
# Fisher Z tranformation for sampled correlations
FisherZ <- function(r){.5*log((1+r)/(1-r))}

# Approximate mean and coviance matrix of Fisher transformed correlations
# having on a IW-distributed covariance matrix. Used in BFcorr and BFcorrUpdate
approxFisherTrans <- function(df,S,samsize=1e4){
  # df is degrees of freedom of inverse-Wishart
  # S is scale matrix of inverse-Wishart
  # samsize is number of draws from IW-distribution to use for approximation
  Sinv <- solve(S)
  drawsPrec <- rWishart(n=samsize,df=df,Sigma=Sinv)
  drawsFC <- sapply(1:samsize,function(s){
    Sigma <- solve(drawsPrec[,,s])
    sds <- sqrt(diag(Sigma))
    CorrMat <- diag(1/sds)%*%Sigma%*%diag(1/sds)
    FisherZ(CorrMat[lower.tri(CorrMat)])
  })
  meanvec <- apply(drawsFC,1,mean)
  covmatrix <- cov(t(drawsFC))
  return(list(mean=meanvec,Sigma=covmatrix))
}

# The function will be called when called BFtest with parameter="correlation"
BFcorr <- function(model,prior="default",constraints="exploratory",priorprob="default"){
  
  #check if 'model' is a fitted mlm-object or list thereof
  if(is.list(model)){
    numpop <- length(model)
    for(pop in 1:numpop){
      if(class(model[[pop]])[1] != "mlm"){stop("Error: 'model' should be a fitted mlm-object with multiple outcome variables or a list of fitted mlm-objects.")}
    }
  }else{ #place the single model in a list object of 1 element
    if(class(model)[1] != "mlm"){stop("Error: 'model' should be a fitted mlm-object with multiple outcome variables or a list of fitted mlm-objects.")}
    numpop <- 1
    model <- list(model)
  }

  # number of outcome varibles
  P <- ncol(model[[1]]$coefficients)
  
  # prior degrees of freedom
  if(!(prior == "default" || prior > 0)){
    stop("Please specify a positive value for the argument 'prior' or leave open and use default")
  }
  if(prior=="default"){
    delta <- 10
  }else{ delta <- prior}
  nu0 <- delta + P - 1
  sd0 <- 1/sqrt(delta + 1) #sd of marginal prior of a correlation
  
  # compute prior mean and covariance matrix of Fisher transformed correlations
  FishApprox <- approxFisherTrans(df=nu0,S=diag(P))
  mean0 <- meanN <- rep(FishApprox$mean,numpop)
  covm0 <- covmN <- kronecker(diag(numpop),FishApprox$Sigma)
  
  # compute posterior mean and covariance matrix of Fisher transformed correlations
  n <- SumSquares <- list()
  for(m in 1:numpop){
    # extract sample size
    n[[m]] <- length(model[[m]]$residuals)
    # compute sums of squares matrix
    SumSquares[[m]] <- cov(model[[m]]$resid)*(n[[m]]-1)
    # posterior mean and covariance matrix of Fisher transformed correlations
    FishApprox <- approxFisherTrans(df=nu0+n[[m]],S=diag(P)+SumSquares[[m]])
    meanN[1:P+(m-1)*P] <- FishApprox$mean
    covmN[1:P+(m-1)*P,1:P+(m-1)*P] <- FishApprox$Sigma
  }
  
  # compute BFs and posterior probs using
  # prior mean en covmatrix mean0 and covm0
  # post mean en covmatrix meanN and covmN
  if(constraints=="exploratory"){
    # H0: corr = 0
    # H1: corr < 0
    # H2: corr < 0
    relfit <- matrix(c(dnorm(0,mean=meanN,sd=sqrt(diag(covmN))),
                       pnorm(0,mean=meanN,sd=sqrt(diag(covmN))),
                       1-pnorm(0,mean=meanN,sd=sqrt(diag(covmN)))),ncol=3)

    relcomp <- matrix(c(dnorm(0,mean=mean0,sd=sqrt(diag(covm0))),
                        pnorm(0,mean=mean0,sd=sqrt(diag(covm0))),
                        1-pnorm(0,mean=mean0,sd=sqrt(diag(covm0)))),ncol=3)
    BFtu <- relfit / relcomp
    PHP <- round(BFtu / apply(BFtu,1,sum),3)
    BFmatrix <- NULL
    
  }else{
    # confirmatory tests based on input constraints

    # for this example the correlations are in the order:
    # corr(DV2,DV1) in group 1 (note that this is the same as corr(DV1,DV2) in group 1)
    # corr(DV3,DV1) in group 1
    # corr(DV3,DV2) in group 1
    # corr(DV2,DV1) in group 2
    # corr(DV3,DV1) in group 2
    # corr(DV3,DV2) in group 2
    
    # Caspar I think we only need a function that generates the labels of these correlations
    names(mean0) <- names(meanN) <- c("r211","r311","r321","r212","r312","r322")
    # I noticed that the bain-function does not work when the covariance matrices have labels
    
    # Example constraints which can be removed when the constraints are read in.
    constraints <- "r211=r311=r321 & r321>r212;
                    r211<(r311,r321) & r321>r212;
                    r211>(r311,r321)"
    
    # compute relative fit
    results <- bain::bain(meanN,constraints,n=10,Sigma=covmN)
    relfit <- results$fit[1:(nrow(results$fit)-1),c(1,3)]
    
    # compute relative complexity
    results <- bain::bain(mean0,constraints,n=10,Sigma=covm0)
    relcomp <- results$fit[1:(nrow(results$fit)-1),c(1,3)]
    
    # get relative fit and complexity of complement hypothesis
    relmeasures <- prob_Hc(mean,covm,relcomp,relfit,constraints)
    relcomp <- relmeasures[[1]]
    relfit <- relmeasures[[2]]
    
    # the BF for the complement hypothesis vs Hu needs to be computed.
    BFtu <- c(apply(relfit / relcomp, 1, prod))
    # Check input of prior probabilies
    if(!(priorprob == "default" || (length(priorprob)==nrow(relfit) && min(priorprob)>0) )){
      stop("'probprob' must be a vector of positive values or set to 'default'.")
    }
    # Change prior probs in case of default setting
    if(priorprob=="default"){priorprobs <- rep(1,length(BFtu))}
    PHP <- round(BFtu*priorprobs / sum(BFtu*priorprobs),3)
    BFmatrix <- matrix(rep(BFtu,length(BFtu)),ncol=length(BFtu))/
      t(matrix(rep(BFtu,length(BFtu)),ncol=length(BFtu)))
  }
  
  return(list(BFtu=BFtu,PHP=PHP,BFmatrix=BFmatrix,relfit=relfit,relcomp=relcomp,
              SumSquares=SumSquares,n=n,constraints=constraints,nu0=nu0,mean0=mean0,
              covm0=covm0,priorprob=priorprob))
}

# The function computes the probability of an unconstrained draw falling in the complement subspace.
prob_Hc <- function(meanvec,covm,relcomp,relfit,constraints){

  numpara <- length(meanvec)
  numhyp <- nrow(relfit)
  relcomp <- relcomp[1:numhyp,]
  which_eq <- relcomp$Fit_eq != 1 
  if(sum(which_eq)==numhyp){ # Then the complement is equivalent to the unconstrained hypothesis.
    relfit <- rbind(relfit,rep(1,2))
    relcomp <- rbind(relcomp,rep(1,2))
    rownames(relfit)[numhyp+1] <- "Hc"
    rownames(relcomp)[numhyp+1] <- "Hc"
  }else{ # So there is at least one hypothesis with only order constraints
    welk <- which(!which_eq)
    if(length(welk)==1){ # There is one hypothesis with only order constraints. Hc is complement of this hypothesis.
      relfit <- rbind(relfit,rep(1,2))
      relfit[numhyp+1,2] <- 1 - relfit[welk,2]
      rownames(relfit)[numhyp+1] <- "Hc"
      relcomp <- rbind(relcomp,rep(1,2))
      relcomp[numhyp+1,2] <- 1 - relcomp[welk,2]
      rownames(relcomp)[numhyp+1] <- "Hc"
    }else{ # So more than one hypothesis with only order constraints
      # First we check whether ther is an overlap between the order constrained spaces.
      
      # Caspar, here we need the RE and RO which are lists of
      # matrices for equality and order constraints under the hypotheses. We can probably do this
      # using the R-code you wrote and a vector of names of the correlations but I don't know
      # how exactly. When running your function I also get an error message saying that he
      # does not know the function "rename_function".
      
      # Maybe (?) create_matrices(varnames=names(meanvec),hyp=constraints)
      
      # This part can be removed when the constraints can be read in.
      constraints <- "r211=r311=r321 & r321>r212;
                      r211<(r311,r321) & r321>r212;
                      r211>(r311,r321)"
      # the coefficient matrices are
      RE1 <- matrix(0,ncol=7,nrow=2)
      RE1[1,c(1,2)] <- c(1,-1)
      RE1[2,c(2,3)] <- c(1,-1)
      RE <- list(RE1,NULL,NULL)
      
      RO1 <- matrix(0,ncol=7,nrow=1)
      RO1[1,c(3,4)] <- c(1,-1)
      RO2 <- matrix(0,ncol=7,nrow=3)
      RO2[1,c(1,2)] <- c(-1,1)
      RO2[2,c(1,3)] <- c(-1,1)
      RO2[3,c(3,4)] <- c(1,-1)
      RO3 <- matrix(0,ncol=7,nrow=2)
      RO3[1,c(1,2)] <- c(1,-1)
      RO3[2,c(1,3)] <- c(1,-1)
      RO <- list(RO1,RO2,RO3)
      #######
      
      draws2 <- 1e4
      priorDraws <- mvtnorm::rmvnorm(draws2,mean=mean0,sigma=covm0)
      #get draws that satisfy the constraints of the separate order constrained hypotheses
      checksOC <- lapply(welk,function(h){
        Rorder <- as.matrix(RO[[h]][,-(1+numpara)])
        if(ncol(Rorder)==1){
          Rorder <- t(Rorder)
        }
        rorder <- as.matrix(RO[[h]][,1+numpara])
        apply(priorDraws%*%t(Rorder) > rep(1,draws2)%*%t(rorder),1,prod)
      })
      checkOCplus <- Reduce("+",checksOC)
      
      if(sum(checkOCplus > 0) < draws2){ #then the joint order constrained hypotheses do not completely cover the parameter space.
        if(sum(checkOCplus>1)==0){ # then order constrained spaces are nonoverlapping
          relfit <- rbind(relfit,rep(1,2))
          relfit[numhyp+1,2] <- 1 - sum(relfit[welk,2])
          rownames(relfit)[numhyp+1] <- "Hc"
          relcomp <- rbind(relcomp,rep(1,2))
          relcomp[numhyp+1,2] <- 1 - sum(relcomp[welk,2])
          rownames(relcomp)[numhyp+1] <- "Hc"
        }else{ #the order constrained subspaces at least partly overlap
          
          # funtion below gives a rough estimate of the posterior probability under Hc
          # a bain type of algorithm would be better of course. but for now this is ok.
          
          relcomp <- rbind(relcomp,rep(1,2))
          relcomp[numhyp+1,2] <- sum(checkOCplus==0)/draws2
          rownames(relcomp)[numhyp+1] <- "Hc"
          
          postDraws <- mvtnorm::rmvnorm(draws2,mean=meanN,sigma=covmN)
          checksOCpost <- lapply(welk,function(h){
            Rorder <- as.matrix(RO[[h]][,-(1+numpara)])
            if(ncol(Rorder)==1){
              Rorder <- t(Rorder)
            }
            rorder <- as.matrix(RO[[h]][,1+numpara])
            apply(postDraws%*%t(Rorder) > rep(1,draws2)%*%t(rorder),1,prod)
          })
          relfit <- rbind(relfit,rep(1,2))
          relfit[numhyp+1,2] <- sum(Reduce("+",checksOCpost) == 0) / draws2
          rownames(relfit)[numhyp+1] <- "Hc"
        }
      }
    }
  }

  return(list(relcomp,relfit))
}

# The update function for BFcorr of new data in fitted object 'model'
BFcorrUpdate <- function(BFcorr1,model){
  
  #Check of dimensions of outcome variables and number of populations match
  if(is.list(model)){
    numpop <- length(model)
    for(pop in 1:numpop){
      if(class(model[[pop]])[1] != "mlm"){stop("Error: 'model' should be a fitted mlm-object with multiple outcome variables or a list of fitted mlm-objects.")}
    }
  }else{ #place the single model in a list object of 1 element
    if(class(model)[1] != "mlm"){stop("Error: 'model' should be a fitted mlm-object with multiple outcome variables or a list of fitted mlm-objects.")}
    numpop <- 1
    model <- list(model)
  }
  if(numpop!=length(BFcorr1$SumSquares)){
    stop("Number of populations of model doesn't match with number in BFcorr object")
  }
  
  # number of outcome varibles
  P <- ncol(model[[1]]$coefficients)
  
  if(P!=ncol(BFcorr1$SumSquares[[1]])){
    stop("Number of dimensions of model doesn't match with number in BFcorr object")
  }
  
  relcomp <- BFcorr1$relcomp
  SumSquares <- BFcorr1$SumSquares
  n <- BFcorr1$n
  nu0 <- BFcorr1$nu0
  mean0 <- meanN <- BFcorr1$mean0
  covm0 <- covmN <- BFcorr1$covm0
  constraints <- BFcorr1$constraints
  priorprob <- BFcorr1$priorprob
  
  # compute posterior mean and covariance matrix of Fisher transformed correlations
  for(m in 1:numpop){
    # extract sample size
    nnew_m <- length(model[[m]]$residuals)
    n[[m]] <- n[[m]] + nnew_m
    # compute sums of squares matrix
    SumSquares[[m]] <- SumSquares[[m]] + cov(model[[m]]$resid)*(nnew_m-1)
    # posterior mean and covariance matrix of Fisher transformed correlations
    FishApprox <- approxFisherTrans(df=nu0+n[[m]],S=diag(P)+SumSquares[[m]])
    meanN[1:P+(m-1)*P] <- FishApprox$mean
    covmN[1:P+(m-1)*P,1:P+(m-1)*P] <- FishApprox$Sigma
  }
  
  # compute BFs and posterior probs using
  # prior mean en covmatrix mean0 and covm0
  # post mean en covmatrix meanN and covmN
  if(constraints=="exploratory"){
    # H0: corr = 0
    # H1: corr < 0
    # H2: corr < 0
    relfit <- matrix(c(dnorm(0,mean=meanN,sd=sqrt(diag(covmN))),
                       pnorm(0,mean=meanN,sd=sqrt(diag(covmN))),
                       1-pnorm(0,mean=meanN,sd=sqrt(diag(covmN)))),ncol=3)
    # relcomp is already available
    BFtu <- relfit / relcomp
    PHP <- round(BFtu / apply(BFtu,1,sum),3)
    BFmatrix <- NULL
    
  }else{
    # compute relative fit
    results <- bain::bain(meanN,constraints,n=10,Sigma=covmN)
    relfit <- results$fit[1:numhyp,c(1,3)]
    # relcomp is already available
    # get relative fit of complement hypothesis
    relmeasures <- prob_Hc(mean,covm,relcomp,relfit,constraints)
    relcomp <- relmeasures[[1]]
    relfit <- relmeasures[[2]]

    # the BF for the complement hypothesis vs Hu needs to be computed.
    BFtu <- c(apply(relfit / relcomp, 1, prod))
    # Change prior probs in case of default setting
    if(priorprob=="default"){priorprobs <- rep(1,length(BFtu))}
    PHP <- round(BFtu*priorprobs / sum(BFtu*priorprobs),3)
    BFmatrix <- matrix(rep(BFtu,length(BFtu)),ncol=length(BFtu))/
      t(matrix(rep(BFtu,length(BFtu)),ncol=length(BFtu)))
    
  }
  
  return(list(BFtu=BFtu,PHP=PHP,BFmatrix=BFmatrix,relfit=relfit,relcomp=relcomp,
              SumSquares=SumSquares,n=n,constraints=constraints,nu0=nu0,mean0=mean0,
              covm0=covm0))
  
}

# Combine two (multivariate) normal distributions (for BFupdate in combination with bain)
# Note. Not used in BFcorr
combineNormals <- function(mean1,covm1,mean2,covm2){
  covm12 <- solve(solve(covm1)+solve(covm2))
  if(nrow(covm12)==1){covm12 <- covm12[1,1]}
  mean12 <- c(covm12%*%(solve(covm1)%*%mean1+solve(covm2)%*%mean2))
  return(list(mean=mean12,covm=covm12))
}



# Example (dummy) correlation test
dt <- as.data.frame(scale(mtcars))
n <- nrow(mtcars)

model1 <- lm(cbind(mpg,cyl,hp) ~ disp + wt, data = dt)
model2 <- lm(cbind(mpg,cyl,hp) ~ carb + gear, data = dt)
model <- list(model1,model2)

BFcorr1 <- BFcorr(model,prior="default",constraints="exploratory",priorprob="default")
BFcorr2 <- BFcorrUpdate(BFcorr1,model) #update with same information new data contained in 'model' (just for the exercise)







