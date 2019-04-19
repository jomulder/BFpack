



# BF test for a univariate normal linear regression type model
BFregression <- function(lm1,constraints="exploratory",priorprob="default"){
  # default BF on location parameters in a univarite normal linear model
  
  K <- length(lm1$coefficients)
  N <- length(lm1$residuals)
  P <- 1 # dimension of outcome variable
  
  Xmat <- model.matrix(lm1)
  yvec <- c(model.matrix(lm1)%*%lm1$coefficients + lm1$residuals)
  dummyX <- rep(F,K)
  for(k in 1:K){ # Check which are dummy variables corresponding to (adjusted) mean parameters
    uniquek <- sort(unique(Xmat[,k]))
    # if(length(uniquek)==2){
    #   if(uniquek[1]==0 && uniquek[2]==1){dummyX[k]<-T} #(adjusted) mean
    # }else{
    #   if(length(uniquek)==1){dummyX[k]<-T} #intercept parameter
    # }
    if(length(uniquek)<=2){dummyX[k]<-T} #group index of intercept
  }
  #number of groups on variations of dummy combinations 
  groupcode <- unique(Xmat[,dummyX])
  rownames(groupcode) <- NULL
  J <- nrow(groupcode)
  if(J==0){ #then no intercept
    J <- 1
    nointercept <- T
  }
  #standardize X
  Xmat[,!dummyX] <- scale(Xmat[,!dummyX])
  # group membership of each observation
  dvec <- unlist(lapply(1:N,function(i){
    which(rowSums(abs(t(matrix(rep(Xmat[i,dummyX],J),ncol=J)) - groupcode))==0)
  }))
  Nj <- c(table(dvec))
  #set minimal fractions for each group
  bj <- ((1+K)/J)/Nj
  # Xj <- lapply(1:J,function(j){
  #   Xmat[dvec==j,]
  #   })
  # yj <- lapply(1:J,function(j){
  #   yvec[dvec==j]
  # })
  #Compute sufficient statistics for all groups
  tXXj <- lapply(1:J,function(j){
    t(Xmat[dvec==j,])%*%Xmat[dvec==j,]
  })
  tXXj_b <- lapply(1:J,function(j){
    t(Xmat[dvec==j,])%*%Xmat[dvec==j,]*bj[j]
  })
  tXyj <- lapply(1:J,function(j){
    t(Xmat[dvec==j,])%*%yvec[dvec==j]
  })
  tXyj_b <- lapply(1:J,function(j){
    t(Xmat[dvec==j,])%*%yvec[dvec==j]*bj[j]
  })
  tyyj <- lapply(1:J,function(j){
    t(yvec[dvec==j])%*%yvec[dvec==j]
  })
  tyyj_b <- lapply(1:J,function(j){
    t(yvec[dvec==j])%*%yvec[dvec==j]*bj[j]
  })
  tXX <- Reduce("+",tXXj)
  tXy <- Reduce("+",tXyj)
  tyy <- Reduce("+",tyyj)
  tXX_b <- Reduce("+",tXXj_b)
  tXy_b <- Reduce("+",tXyj_b)
  tyy_b <- Reduce("+",tyyj_b)
  betaHat <- solve(tXX)%*%tXy           # same as lm1$coefficients
  s2 <- tyy - t(tXy)%*%solve(tXX)%*%tXy # same as sum((lm1$residuals)**2)
  # sufficient statistics based on fraction of the data
  betaHat_b <- solve(tXX_b)%*%tXy_b
  s2_b <- tyy_b - t(tXy_b)%*%solve(tXX_b)%*%tXy_b
  
  if(constraints=="exploratory"){
    
    # prior hyperparameters
    scale0 <- c(s2_b)*solve(tXX_b)
    df0 <- 1
    mean0 <- as.matrix(rep(0,K))
    # posterior hyperparameters
    scaleN <- c(s2)*solve(tXX)/(N-K)
    dfN <- N-K
    meanN <- betaHat
    
    # Hypotheses for exploraotyr test
    # H0: beta = 0
    # H1: beta < 0
    # H2: beta < 0
    relfit <- t(matrix(unlist(lapply(1:K,function(k){
        c(dt((0-meanN[k,1])/sqrt(scaleN[k,k]),df=dfN)/sqrt(scaleN[k,k]),
                           pt((0-meanN[k,1])/sqrt(scaleN[k,k]),df=dfN,lower.tail = TRUE),
                           pt((0-meanN[k,1])/sqrt(scaleN[k,k]),df=dfN,lower.tail = FALSE))
      })),nrow=3))
    relcomp <- t(matrix(unlist(lapply(1:K,function(k){
        c(dt((0-mean0[k,1])/sqrt(scale0[k,k]),df=df0)/sqrt(scale0[k,k]),
          pt((0-mean0[k,1])/sqrt(scale0[k,k]),df=df0,lower.tail = TRUE),
          pt((0-mean0[k,1])/sqrt(scale0[k,k]),df=df0,lower.tail = FALSE))
      })),nrow=3))
    
    BFtu <- relfit / relcomp
    PHP <- round(BFtu / apply(BFtu,1,sum),3)
    BFmatrix <- NULL

  }else{
  
    #read constraints. FOR NOW AN EXAMPLE. Use Caspar's function here.
    RrE1 <- matrix(0,ncol=7,nrow=2)
    RrE1[1,c(4,5)] <- c(1,-1)
    RrE1[2,c(4,6)] <- c(1,-1)
    RrE <- list(RrE1,RrE1,NULL)
    
    RrO1 <- matrix(0,ncol=7,nrow=1)
    RrO1[1,c(3,4)] <- c(1,-1)
    RrO2 <- matrix(0,ncol=7,nrow=3)
    RrO2[1,c(1,2)] <- c(-1,1)
    RrO2[2,c(1,3)] <- c(-1,1)
    RrO2[3,c(3,4)] <- c(1,-1)
    RrO3 <- matrix(0,ncol=7,nrow=2)
    RrO3[1,c(1,2)] <- c(1,-1)
    RrO3[2,c(1,3)] <- c(1,-1)
    RrO <- list(RrO1,NULL,RrO2)
    #######
    
    RrStack <- rbind(do.call(rbind,RrE),do.call(rbind,RrO))
    RStack <- RrStack[,-(K+1)]
    rStack <- RrStack[,(K+1)]
    
    # check if a common boundary exists for prior location under all constrained hypotheses
    if(nrow(RrStack) > 1){
      rref_ei <- pracma::rref(RrStack)
      nonzero <- RrStack[,K+1]!=0
      if(max(nonzero)>0){
        row1 <- max(which(nonzero==T))
        if(sum(abs(RrStack[row1,1:K]))==0){
          stop("No common boundary point for prior location. Conflicting constraints.")
        }
      }
    }
    # prior hyperparameters
    scale0 <- c(s2_b)*solve(tXX_b)
    df0 <- 1
    mean0 <- MASS::ginv(RStack)%*%rStack
    # posterior hyperparameters
    scaleN <- c(s2)*solve(tXX)/(N-K)
    dfN <- N-K
    meanN <- betaHat
    numhyp <- length(RrO)
    
    relcomp <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Student_measures(mean0,scale0,df0,RrE[[h]],RrO[[h]])
    })),nrow=2))
    
    relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Student_measures(meanN,scaleN,dfN,RrE[[h]],RrO[[h]])
    })),nrow=2))
    
    # ######## ### ######## # ##### # ###
    # Use complement hypothesis function from correlation_testing
    # Alter this function such that the prob is computed given a N(m,Sigma)
    # and a t(m,S,df) instead of giving both relcomp and relfit in one run.
    # ###### ########## # ## ###### #####
    
    relfit <- Student_prob_Hc(meanN,scaleN,dfN,relfit,constraints)
    relcomp <- Student_prob_Hc(mean0,scale0,df0,relcomp,constraints)
    
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
              Nj=Nj,bj=bj,tXXj=tXXj,tXyj=tXyj,tyyj=tyyj,constraints=constraints,
              priorprob=priorprob))
}

# compute relative meausures (fit or complexity) under a multivariate Student t distribution
Student_measures <- function(mean1,scale1,df1,RrE1,RrO1){
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
    scaleE <- RE1%*%scale1%*%t(RE1)
    relE <- dmvt(rE1,delta=c(meanE),sigma=scaleE,df=df1,log=FALSE)
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
      scaleO <- RO1%*%scale1%*%t(RO1)
      relO <- ifelse(nrow(scaleO)==1,
                pt((rO1-meanO)/sqrt(scaleO[1,1]),df=df1,lower.tail=FALSE), #univariate
                pmvt(lower=rO1,upper=Inf,delta=meanO,sigma=scaleO,df=df1,type="shifted")) #multivariate
    }else{ #no linear transformation can be used; pmvt cannot be used. Use bain with a multivariate normal approximation
      #compute covariance matrix for multivariate normal distribution
      Sigma1 <- ifelse(df1>2,
                  df1/(df1-2)*scale1,
                  scale1) #for prior with df1==1, probability independent of common factor of scale1
      # bain1 <- bain::bain(mean1,Sigma1=Sigma1,RrE1,RrO1,n=10) # choice of n does not matter
      # extract posterior probability (Fit_eq) from bain-object)
      stop("REMINDER. This situation should still be implemented.")
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
    
    #a)Transformation matrix
    D <- diag(K) - t(RE1) %*% solve(RE1 %*% t(RE1)) %*% RE1
    D2 <- unique(round(D, 5))
    D2 <- D2[as.logical(rowSums(D2 != 0)),]
    Tm <- rbind(RE1, D2)
    
    #b)
    Tmean1 <- Tm %*% mean1
    Tscale1 <- Tm %*% scale1 %*% t(Tm)
    
    # relative meausure for equalities
    relE <- mvtnorm::dmvt(x = t(rE1), delta = Tmean1[1:qE1], sigma = matrix(Tscale1[1:qE1, 1:qE1], ncol = qE1), df = N - K, log = FALSE)
    
    # transform order constraints
    RO1tilde <- RO1 %*% MASS::ginv(D2)
    rO1tilde <- rO1 - RO1 %*% MASS::ginv(RE1) %*% rE1
    
    # Partitioning equality part and order part
    Tmean1E <- Tmean1[1:qE1]
    Tmean1O <- Tmean1[(qE1 + 1):K]
    
    Tscale1EE <- Tscale1[1:qE1, 1:qE1]
    Tscale1OE <- Tscale1[(qE1 + 1):K, 1:qE1]
    Tscale1OO <- Tscale1[(qE1 + 1):K, (qE1 + 1):K]
    
    #conditional location and scale matrix
    Tmean1OgE <- Tmean1O + Tscale1OE %*% solve(Tscale1EE) %*% matrix(rE1 - Tmean1E)
    Tscale1OgE <- as.vector((N - K + (t(matrix(rE1 - Tmean1E)) %*% solve(Tscale1EE) %*% matrix(rE1 - Tmean1E))) /
                              (N - K + qE1)) * (Tscale1OO - Tscale1OE %*% solve(Tscale1EE) %*% t(Tscale1OE))
    
    if(Matrix::rankMatrix(RO1tilde)[[1]] == nrow(RO1tilde)){
      rO1tilde <- as.vector(rO1tilde)
      
      delta_trans <- as.vector(RO1tilde %*% Tmean1OgE)
      scale1_trans <- RO1tilde %*% Tscale1OgE %*% t(RO1tilde)
      
      if(nrow(scale1_trans) == 1){ # univariate
        relO <- pt((rO1tilde - delta_trans) / sqrt(scale1_trans), df = N-K+qE1, lower.tail = FALSE)[1]
      } else { # multivariate
        relO <- mvtnorm::pmvt(lower = rO1tilde, upper = Inf, delta = delta_trans, sigma = scale1_trans, df = N - K + qE1, type = "shifted")[1]
      }
      
    }else{ #use bain for the computation of the probability
      #compute covariance matrix for multivariate normal distribution
      Sigma1 <- ifelse(df1>2,
                       df1/(df1-2)*Tscale1OgE,
                       Tscale1OgE) #for prior with df1==1, probability independent of common factor of scale1
      # bain1 <- bain::bain(mean1,Sigma1=Sigma1,RrE1=NULL,RrO1=RO1tilde,n=10) # choice of n does not matter
      # extract posterior probability (Fit_eq) from bain-object)
      stop("REMINDER. This situation should still be implemented.")
    }
  }
  
  return(c(relE,relO))
}

# The function computes the probability of an unconstrained draw falling in the complement subspace.
Student_prob_Hc <- function(mean1,scale1,df1,relmeas,constraints){
  
  numpara <- length(mean1)
  numhyp <- nrow(relmeas)
  relmeas <- relmeas[1:numhyp,]
  which_eq <- relmeas[,1] != 1 
  if(sum(which_eq)==numhyp){ # Then the complement is equivalent to the unconstrained hypothesis.
    relmeas <- rbind(relmeas,rep(1,2))
    rownames(relfit)[numhyp+1] <- "Hc"
  }else{ # So there is at least one hypothesis with only order constraints
    welk <- which(!which_eq)
    if(length(welk)==1){ # There is one hypothesis with only order constraints. Hc is complement of this hypothesis.
      relmeas <- rbind(relmeas,rep(1,2))
      relmeas[numhyp+1,2] <- 1 - relmeas[welk,2]
      rownames(relmeas)[numhyp+1] <- "Hc"
    }else{ # So more than one hypothesis with only order constraints
      # First we check whether ther is an overlap between the order constrained spaces.
      
      # Caspar, here we need the RE and RO which are lists of
      # matrices for equality and order constraints under the hypotheses. We can probably do this
      # using the R-code you wrote and a vector of names of the correlations but I don't know
      # how exactly. When running your function I also get an error message saying that he
      # does not know the function "rename_function".
      
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
          
          randomDraws <- mvtnorm::rmvt(draws2,delta=mean1,sigma=scale1,df=df1)
          checksOC <- lapply(welk,function(h){
            Rorder <- as.matrix(RrO[[h]][,-(1+numpara)])
            if(ncol(Rorder)==1){
              Rorder <- t(Rorder)
            }
            rorder <- as.matrix(RrO[[h]][,1+numpara])
            apply(randomDraws%*%t(Rorder) > rep(1,draws2)%*%t(rorder),1,prod)
          })
          relmeas <- rbind(relmeas,rep(1,2))
          relmeas[numhyp+1,2] <- sum(Reduce("+",checksOC) == 0) / draws2
          rownames(relmeas)[numhyp+1] <- "Hc"
        }
      }
    }
  }
  
  return(relmeas)
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
    relE <- dmvnorm(rE1,mean=c(meanE),sigma=SigmaE,log=FALSE)
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
      relO <- pmvt(lower=rO1,upper=Inf,mean=meanO,sigma=SigmaO)
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
    
    #a)Transformation matrix
    D <- diag(K) - t(RE1) %*% solve(RE1 %*% t(RE1)) %*% RE1
    D2 <- unique(round(D, 5))
    D2 <- D2[as.logical(rowSums(D2 != 0)),]
    Tm <- rbind(RE1, D2)
    
    #b)
    Tmean1 <- Tm %*% mean1
    TSigma1 <- Tm %*% Sigma1 %*% t(Tm)
    
    # relative meausure for equalities
    relE <- mvtnorm::dmvnorm(x=rE1,mean=Tmean1[1:qE1],sigma=matrix(TSigma1[1:qE1,1:qE1],ncol=qE1),log=FALSE)
    
    # transform order constraints
    RO1tilde <- RO1 %*% MASS::ginv(D2)
    rO1tilde <- rO1 - RO1 %*% MASS::ginv(RE1) %*% rE1
    
    # Partitioning equality part and order part
    Tmean1E <- Tmean1[1:qE1]
    Tmean1O <- Tmean1[(qE1+1):K]
    
    TSigma1EE <- TSigma1[1:qE1,1:qE1]
    TSigma1OE <- TSigma1[(qE1+1):K,1:qE1]
    TSigma1OO <- TSigma1[(qE1+1):K,(qE1+1):K]
    
    #conditional location and covariance matrix
    Tmean1OgE <- Tmean1O + TSigma1OE %*% solve(TSigma1EE) %*% (rE1-Tmean1E)
    TSigma1OgE <- TSigma1OO - TSigma1OE %*% solve(TSigma1EE) %*% t(TSigma1OE)
    
    if(Matrix::rankMatrix(RO1tilde)[[1]] == nrow(RO1tilde)){
      rO1tilde <- as.vector(rO1tilde)
      mean_trans <- as.vector(RO1tilde %*% Tmean1OgE)
      Sigma1_trans <- RO1tilde %*% TSigma1OgE %*% t(RO1tilde)
      
      relO <- mvtnorm::pmvnorm(lower=rO1tilde,upper=Inf,mean=mean_trans,sigma=Sigma1_trans,log=FALSE)
      
    }else{ #use bain for the computation of the probability
      
      # bain1 <- bain::bain(mean1,Sigma1=Sigma1,RrE1=NULL,RrO1=RO1tilde,n=10) # choice of n does not matter
      # extract posterior probability (Fit_eq) from bain-object)
      stop("REMINDER. This case should still be implemented.")
    }
  }
  
  return(c(relE,relO))
}



mtcars0 <- mtcars
mtcars0$vs[1:6] <- 2
mtcars0$vs <- as.factor(mtcars0$vs)
mtcars0$am <- as.factor(mtcars0$am)
summary(mtcars0)

lm1 <- lm(wt ~ -1 + disp + vs + hp + drat, mtcars0)
model.matrix(lm1)
summary(lm1)

BFregression(lm1,constraints="exploratory")
























