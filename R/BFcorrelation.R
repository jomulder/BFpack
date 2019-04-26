
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
BFcorr <- function(model,prior=NULL,constraints="exploratory",priorprob="default"){

  #check if 'model' is a fitted mlm-object or list thereof
  if(class(model)[1] == "mlm"){
    numgroup <- 1
    model <- list(model)
  }else{ if(is.list(model)){
    numgroup <- length(model)
  }else{
    stop("Error: 'model' should be a fitted mlm-object with multiple outcome variables or a list thereof.")
  }
  }

  # number of outcome varibles
  P <- ncol(model[[1]]$coefficients)
  numcorrgroup <- P*(P-1)/2

  # prior degrees of freedom
  if(is.null(prior)){
    delta <- 10
  }else if(is.numeric(prior)){
    if(prior>0){
      delta <- prior
    }else stop("Please specify a positive value for the argument 'prior' or use default")
  }else stop("Please specify a positive value for the argument 'prior' or use default")
  nu0 <- delta + P - 1
  sd0 <- 1/sqrt(delta + 1) #sd of marginal prior of a correlation

  # compute prior mean and covariance matrix of Fisher transformed correlations
  FishApprox <- approxFisherTrans(df=nu0,S=diag(P))
  mean0 <- meanN <- rep(0,numcorrgroup*numgroup)
  covm0 <- covmN <- kronecker(diag(numgroup),FishApprox$Sigma)

  # compute posterior mean and covariance matrix of Fisher transformed correlations
  n <- SumSquares <- list()
  for(m in 1:numgroup){
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

    corrmat <- diag(P)
    row.names(corrmat) <- colnames(corrmat) <- row.names(SumSquares[[m]])
    corr_names <- names(get_estimates(corrmat)$estimate)
    matrix_names <- matrix(corr_names,nrow=3)
    # equal correlations are at the opposite side of the vector
    corr_names <- c(matrix_names[lower.tri(matrix_names)],
      t(matrix_names)[lower.tri(matrix_names)])
    if(numgroup>1){
      matcorrpop <- matrix(0,nrow=length(corr_names),ncol=numgroup)
      for(c in 1:length(corr_names)){
        matcorrpop[c,] <- unlist(lapply(1:numgroup,function(pop){
          paste0(corr_names[c],"_gr",as.character(pop)) #or "_in_gr"
        }))
      }
      corr_names <- c(matcorrpop)
    }

    parse_hyp <- parse_hypothesis(corr_names,constraints)
    select1 <- rep(1:numcorrgroup,numgroup) + rep((0:(numgroup-1))*2*numcorrgroup,each=numcorrgroup)
    select2 <- rep(numcorrgroup+1:numcorrgroup,numgroup) + rep((0:(numgroup-1))*2*numcorrgroup,each=numcorrgroup)
    #combine equivalent correlations, e.g., cor(Y1,Y2)=corr(Y2,Y1).
    parse_hyp$hyp_mat <-
      cbind(parse_hyp$hyp_mat[,select1] + parse_hyp$hyp_mat[,select2],parse_hyp$hyp_mat[,numcorrgroup*2*numgroup+1])
    #create coefficient with equality and order constraints
    RrList <- make_RrList(parse_hyp)
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]

    numhyp <- length(RrE)
    relcomp <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Gaussian_measures(mean0,covm0,RrE[[h]],RrO[[h]])
    })),nrow=2))
    relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Gaussian_measures(meanN,covmN,RrE[[h]],RrO[[h]])
    })),nrow=2))

    # get relative fit and complexity of complement hypothesis
    relcomp <- Gaussian_prob_Hc(mean0,covm0,relcomp,RrO)
    relfit <- Gaussian_prob_Hc(meanN,covmN,relfit,RrO)

    Hnames <- c(unlist(lapply(1:numhyp,function(h){paste0("H",as.character(h))})),"Hc")
    row.names(relcomp) <- Hnames
    row.names(relfit) <- Hnames

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
    row.names(BFmatrix) <- Hnames
    colnames(BFmatrix) <- Hnames
  }

  return(list(BFtu=BFtu,PHP=PHP,BFmatrix=BFmatrix,relfit=relfit,relcomp=relcomp,
              SumSquares=SumSquares,n=n,constraints=constraints,nu0=nu0,mean0=mean0,
              covm0=covm0,priorprob=priorprob))
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

# compute relative meausures (fit or complexity) under a multivariate Student t distribution
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
      TSigma1OgE <- TSigma1OO - TSigma1OE %*% solve(TSigma1EE) %*% t(TSigma1OE)

      relO <- mvtnorm::pmvnorm(lower=rO1,upper=Inf,mean=c(Tmean1OgE),sigma=TSigma1OgE)[1]

    }else{ #use bain for the computation of the probability

      # bain1 <- bain::bain(mean1,Sigma1=Sigma1,RrE1=NULL,RrO1=RO1tilde,n=10) # choice of n does not matter
      # extract posterior probability (Fit_eq) from bain-object)
      stop("REMINDER. This case should still be implemented.")
    }
  }

  return(c(relE,relO))
}

# The function computes the probability of an unconstrained draw falling in the complement subspace.
Gaussian_prob_Hc <- function(mean1,Sigma1,relmeas,constraints,RrO){

  numpara <- length(mean1)
  numhyp <- nrow(relmeas)
  relmeas <- relmeas[1:numhyp,]
  which_eq <- relmeas[,1] != 1
  if(sum(which_eq)==numhyp){ # Then the complement is equivalent to the unconstrained hypothesis.
    relmeas <- rbind(relmeas,rep(1,2))
    rownames(relmeas)[relmeas+1] <- "Hc"
  }else{ # So there is at least one hypothesis with only order constraints
    welk <- which(!which_eq)
    if(length(welk)==1){ # There is one hypothesis with only order constraints. Hc is complement of this hypothesis.
      relmeas <- rbind(relmeas,rep(1,2))
      relmeas[numhyp+1,2] <- 1 - relmeas[welk,2]
      rownames(relmeas)[numhyp+1] <- "Hc"
    }else{ # So more than one hypothesis with only order constraints
      # First we check whether ther is an overlap between the order constrained spaces.
#
#       # Caspar, here we need the RE and RO which are lists of
#       # matrices for equality and order constraints under the hypotheses. We can probably do this
#       # using the R-code you wrote and a vector of names of the correlations but I don't know
#       # how exactly. When running your function I also get an error message saying that he
#       # does not know the function "rename_function".
#
#       # Maybe (?) create_matrices(varnames=names(meanvec),hyp=constraints)
#
#       # This part can be removed when the constraints can be read in.
#       constraints <- "r211=r311=r321 & r321>r212;
#       r211<(r311,r321) & r321>r212;
#       r211>(r311,r321)"
#       # the coefficient matrices are
#       RE1 <- matrix(0,ncol=7,nrow=2)
#       RE1[1,c(1,2)] <- c(1,-1)
#       RE1[2,c(2,3)] <- c(1,-1)
#       RE <- list(RE1,NULL,NULL)
#
#       RO1 <- matrix(0,ncol=7,nrow=1)
#       RO1[1,c(3,4)] <- c(1,-1)
#       RO2 <- matrix(0,ncol=7,nrow=3)
#       RO2[1,c(1,2)] <- c(-1,1)
#       RO2[2,c(1,3)] <- c(-1,1)
#       RO2[3,c(3,4)] <- c(1,-1)
#       RO3 <- matrix(0,ncol=7,nrow=2)
#       RO3[1,c(1,2)] <- c(1,-1)
#       RO3[2,c(1,3)] <- c(1,-1)
#       RO <- list(RO1,RO2,RO3)
#       #######

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
          rownames(relmeas)[numhyp+1] <- "Hc"
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
          rownames(relmeas)[numhyp+1] <- "Hc"
        }
      }
    }
  }

  return(relmeas)
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

constraints <- "cyl_with_mpg_gr1 > hp_with_mpg_gr1 > hp_with_cyl_gr1; cyl_with_mpg_gr1 < hp_with_mpg_gr1 = cyl_with_hp_gr1 ; cyl_with_mpg_gr2 = hp_with_mpg_gr2 = cyl_with_hp_gr1"
constraints <- "cyl_with_mpg > hp_with_mpg > hp_with_cyl; cyl_with_mpg < hp_with_mpg = cyl_with_hp ; cyl_with_mpg = hp_with_mpg = cyl_with_hp"

BFcorr1 <- BFcorr(model,prior=NULL,constraints="exploratory",priorprob="default")
BFcorr2 <- BFcorrUpdate(BFcorr1,model) #update with same information new data contained in 'model' (just for the exercise)







