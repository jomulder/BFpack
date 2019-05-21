
# Compute relative meausures for constraints on correlation based
# on a joint uniform distribution
jointuniform_measures <- function(P,numcorrgroup,numG,RrE1,RrO1,draws){
  
  relE <- relO <- 1
  # If the correlation
  teldummy <- 0
  if (numcorr==1){ #use analytic expression in case of a single correlation
    # RrE1=RrO1=matrix(c(-1,1,-.5,.3),nrow=2)
    if(!is.null(RrE1)){ # only inequality constraint(s). Copmute prior proportion.
      relE <- .5
    }else if(!is.null(RrO1)){
      if(nrow(RrO1)==1){
        relO <- .5*ifelse(RrO1[1,1]>0,1-RrO1[1,2]/RrO1[1,1],RrO1[1,2]/RrO1[1,1]+1)
      }else if(nrow(RrO1)==2){
        which_pos <- which(RrO1[,1]>0)
        which_neg <- (1:2)[-which_pos]
        min1 <- RrO1[which_pos,2]/RrO1[which_pos,1]
        max1 <- RrO1[which_neg,2]/RrO1[which_neg,1]
        relO <- (max1-min1)/2
      }else{
        stop("The constraints on the correlations seems to be conflicting under one hypothesis.")
      }
    }
  }else{ #use numerical estimate in case of a multiple correlations
    
    # get random draws from joint uniform distribution
    drawsJU <- matrix(0,nrow=samsize,ncol=numcorr)
    drawsJU[,1:numcorrgroup] <- .Fortran("draw_ju",P=as.integer(P),drawscorr=drawsJU[,1:numcorrgroup],
                                         samsize=as.integer(samsize),numcorrgroup=as.integer(numcorrgroup),
                                         seed=as.integer(seed))$drawscorr
    # REMINDER: Check if Fisher transformed draws are needed or nontransformed draws.

    #These draws can be used for the correlation matrices of the independent groups
    numcorr <- numcorrgroup * numG
    teldummy <- numcorrgroup
    if(numG>1){
      for(i in 2:numG){
        drawsJU[1:(samsize-i+1),(teldummy+1):(teldummy+numcorrgroup)] = drawsJU[i:samsize,1:numcorrgroup]
        drawsJU[(samsize-i+2):samsize,(teldummy+1):(teldummy+numcorrgroup)] = drawsJU[1:(i-1),1:numcorrgroup]
        teldummy <- teldummy + numcorrgroup
      }
    }
    
    if(!is.null(RrE1) && is.null(RrO1)){ #only equality constraints
      
      # RrE1 <- matrix(0,ncol=7,nrow=2)
      # RrE1[1,c(1,4)] <- c(1,-1)
      # RrE1[2,1+c(1,4)] <- c(1,-1)
      
      if(nrow(RrE1)==1){
        RE1 <- t(RrE1[1,-numcorr-1])
        rE1 <- RrE1[1,numcorr+1]
      }else{
        RE1 <- RrE1[,-numcorr-1]
        rE1 <- RrE1[,numcorr+1]
      }
      drawsTrans = drawsJU%*%t(RE1)
      relE <- .Fortran("compute_rcEt",numE=as.integer(nrow(RrE1)),drawsIn=drawsTrans,wIn=rE1,delta=.1,
                       rcEt=relE,samsize=as.integer(samsize))$rcEt
    }else if(is.null(RrE1) && !is.null(RrO1)){ #only order constraints.
      # Use normal approximation.
      relO <- Gaussian_measures(mean1=mean(drawsJU),Sigma1=cov(drawsJU),RrE1=NULL,RrO1)[2]
    }else if(!is.null(RrE1) && !is.null(RrO1)){
      # RrO1 <- matrix(0,ncol=7,nrow=1)
      # RrO1[1,2+c(1,4)] <- c(1,-1)
      if(nrow(RrE1)==1){
        RE1 <- t(RrE1[1,-numcorr-1])
        rE1 <- RrE1[1,numcorr+1]
      }else{
        RE1 <- RrE1[,-numcorr-1]
        rE1 <- RrE1[,numcorr+1]
      }
      if(nrow(RrO1)==1){
        RO1 <- t(RrO1[1,-numcorr-1])
        rO1 <- RrO1[1,numcorr+1]
      }else{
        RO1 <- RrO1[,-numcorr-1]
        rO1 <- RrO1[,numcorr+1]
      }
      if(Matrix::rankMatrix(RrO1)[[1]] == nrow(RrO1)){
        R1 <- rbind(RE1,RO1)
        r1 <- c(rE1,rO1)
        drawsTrans = drawsJU%*%t(R1)
        
        meanOE <- rep(0,nrow(RrO1))
        covOE <- matrix(0,nrow=nrow(RrO1),ncol=nrow(RrO1))
        analysisE <- .Fortran("compute_rcEt2",numE=as.integer(nrow(RrE1)),drawsIn=drawsTrans,wIn=rE1,delta=.1,
                 rcEt=relE,meanOut=meanOE,covmOut=covOE,samsize=as.integer(samsize),numcorr=as.integer(nrow(R1)))
        relE <- analysisE$rcEt
        meanOE <- analysisE$meanOut
        covOE <- analysisE$covmOut
        # Use normal approximation for the conditional probability
        RrOE1 <- cbind(diag(nrow(RrO1)),rO1)
        relO <- Gaussian_measures(mean1=meanOE,Sigma1=covOE,RrE1=NULL,RrO1=RrOE1)[2]
      }else { # transformation with RO1 is not possible
        
        #make one-to-one transformation matrix Tm
        D <- diag(numcorr) - t(RE1) %*% solve(RE1 %*% t(RE1)) %*% RE1
        D2 <- unique(round(D, 5))
        D2 <- D2[as.logical(rowSums(D2 != 0)),]
        Tm <- rbind(RE1, D2)
        
        drawsTrans = drawsJU%*%t(Tm)
        
        meanOE <- rep(0,numcorr-nrow(RE1))
        covOE <- matrix(0,nrow=numcorr-nrow(RE1),ncol=numcorr-nrow(RE1))
        analysisE <- .Fortran("compute_rcEt2",numE=as.integer(nrow(RrE1)),drawsIn=drawsTrans,wIn=rE1,delta=.1,
                          rcEt=relE,meanOut=meanOE,covmOut=covOE,samsize=as.integer(samsize),numcorr=as.integer(numcorr))
        relE <- analysisE$rcEt
        meanOE <- analysisE$meanOut
        covOE <- analysisE$covmOut
        
        ROE1 <- RO1 %*% MASS::ginv(D2)
        rOE1 <- rO1 - RO1 %*% MASS::ginv(RE1) %*% rE1
        
        # Use normal approximation for the conditional probability
        relO <- Gaussian_measures(mean1=meanOE,Sigma1=covOE,RrE1=NULL,RrO1=cbind(ROE1,rOE1))[2]
      }
    }
  }
  return(c(relE,relO))
}





