


cor_test_new <- function(..., formula = NULL, iter = 5e3, burnin = 3e3){

  Y_groups <- list(...)
  numG <- length(Y_groups)
  P <- ncol(Y_groups[[1]])
  ordi <- numcats <- matrix(0,nrow=numG,ncol=P)

  #check class of variables in input data
  for(gg in 1:numG){
    for(pp in 1:P){
      if(class(Y_groups[[gg]][,pp])[1] != "numeric"){
        if(class(Y_groups[[gg]][,pp])[1] == "ordered"){
          levels(Y_groups[[gg]][,pp]) <- 1:length(levels(Y_groups[[gg]][,pp]))
          Y_groups[[gg]][,pp] <- as.numeric(Y_groups[[gg]][,pp])
          ordi[gg,pp] <- 1
          numcats[gg,pp] <- max(Y_groups[[gg]][,pp])
        }else{
          if(class(Y_groups[[gg]][,pp])[1] == "factor"){
            if(length(levels(Y_groups[[gg]][,pp]))==2){
              levels(Y_groups[[gg]][,pp]) <- 1:length(levels(Y_groups[[gg]][,pp]))
              Y_groups[[gg]][,pp] <- as.numeric(Y_groups[[gg]][,pp])
              ordi[gg,pp] <- 1
              numcats[gg,pp] <- 2
            }else{
              stop("Outcome variables should be either of class 'numeric', 'ordered', or a 2-level 'factor'.")
            }
          }else{
            stop("Outcome variables should be either of class 'numeric', 'ordered', or a 2-level 'factor'.")
          }
        }
      }
    }
  }

  if(is.null(formula)){
    formula <- ~ 1
  }

  model_matrices <- lapply(seq_len(numG) , function(x) {
    model.matrix(formula, Y_groups[[x]])
  })

  correlate <- remove_predictors_helper(Y_groups = Y_groups, formula)

  YXlist <- lapply(1:length(model_matrices),function(g){
    list(as.matrix(correlate[[g]]),as.matrix(model_matrices[[g]]))
  })

  K <- ncol(YXlist[[1]][[2]])
  numcorr <- numG*P*(P-1)/2
  ngroups <- unlist(lapply(1:numG,function(g){nrow(YXlist[[g]][[1]])}))
  Ntot <- max(ngroups)
  Ygroups <- array(0,dim=c(numG,Ntot,P))
  Xgroups <- array(0,dim=c(numG,Ntot,K))
  XtXi <- array(0,dim=c(numG,K,K))
  BHat <- array(0,dim=c(numG,K,P))
  sdHat <- matrix(0,nrow=numG,ncol=P)
  CHat <- array(0,dim=c(numG,P,P))
  SumSq <- array(0,dim=c(numG,P,P))
  SumSqInv <- array(0,dim=c(numG,P,P))
  sdsd <- matrix(0,nrow=numG,ncol=P)

  for(g in 1:numG){
    Y_g <- YXlist[[g]][[1]]
    for(p in 1:P){
      if(ordi[g,p]==0){
        Y_g[,p] <- c(scale(Y_g[,p]))
      }
    }
    X_g <- YXlist[[g]][[2]]
    Ygroups[g,1:ngroups[g],] <- Y_g
    #standardize data to get a more stable sampler for the correlations.
    tableX <- apply(X_g,2,table)
    catX <- unlist(lapply(1:length(tableX),function(xcol){
      length(tableX[[xcol]])
    }))
    if(sum(catX>1)){
      X_g[1:ngroups[g],which(catX>1)] <- apply(as.matrix(X_g[1:ngroups[g],which(catX>1)]),2,scale)
    }
    Xgroups[g,1:ngroups[g],] <- X_g
    XtXi[g,,] <- solve(t(X_g)%*%X_g)
    BHat[g,,] <- XtXi[g,,]%*%t(X_g)%*%Y_g
    SumSq[g,,] <- t(Y_g - X_g%*%BHat[g,,])%*%(Y_g - X_g%*%BHat[g,,])
    SumSqInv[g,,] <- solve(SumSq[g,,])
    Sigma_g <- SumSq[g,,]/ngroups[g]
    sdHat[g,] <- sqrt(diag(Sigma_g))
    CHat[g,,] <- diag(1/sdHat[g,])%*%Sigma_g%*%diag(1/sdHat[g,])
    #get rough estimate of posterior sd of the standard deviations (used for random walk sd)
    drawsSigma_g <- rWishart(1e2,df=ngroups[g],Sigma=SumSqInv[g,,])
    sdsd[g,] <- unlist(lapply(1:P,function(p){
      sd(sqrt(drawsSigma_g[p,p,]))
    }))
  }
  samsize0 <- iter
  gLiuSab <- array(0,dim=c(samsize0,numG,P))

  # call Fortran subroutine for Gibbs sampling using noninformative improper priors
  # for regression coefficients, Jeffreys priors for standard deviations, and a proper
  # joint uniform prior for the correlation matrices.

  res <- .Fortran("estimate_bct_ordinal",
                  postZmean=matrix(0,numcorr,1),
                  postZcov=matrix(0,numcorr,numcorr),
                  P=as.integer(P),
                  numcorr=as.integer(numcorr),
                  K=as.integer(K),
                  numG=as.integer(numG),
                  BHat=round(BHat,3),
                  sdHat=round(sdHat,3),
                  CHat=round(CHat,3),
                  XtXi=XtXi,
                  samsize0=as.integer(samsize0),
                  burnin=as.integer(burnin),
                  Ntot=as.integer(Ntot),
                  Xgroups=Xgroups,
                  Ygroups=Ygroups,
                  C_quantiles=array(0,dim=c(numG,P,P,3)),
                  sigma_quantiles=array(0,dim=c(numG,P,3)),
                  B_quantiles=array(0,dim=c(numG,K,P,3)),
                  BDrawsStore=array(0,dim=c(samsize0,numG,K,P)),
                  sigmaDrawsStore=array(0,dim=c(samsize0,numG,P)),
                  CDrawsStore=array(0,dim=c(samsize0,numG,P,P)),
                  sdMH=sdsd,
                  ordinal_in=ordi,
                  Cat_in=numcats,
                  maxCat=as.integer(max(numcats)),
                  gLiuSab=gLiuSab,
                  seed=as.integer( sample.int(1e6,1) ))

  varnames <- lapply(1:numG,function(g){
    names(correlate[[g]])
  })
  corrnames <- lapply(1:numG,function(g){
    matrix(unlist(lapply(1:P,function(p2){
      unlist(lapply(1:P,function(p1){
        if(numG==1){
          paste0(varnames[[g]][p1],"_with_",varnames[[g]][p2])
        }else{
          paste0(varnames[[g]][p1],"_with_",varnames[[g]][p2],"_in_g",as.character(g))
        }
      }))
    })),nrow=P)
  })

  FmeansCovCorr <- lapply(1:numG,function(g){
    Fdraws_g <- FisherZ(t(matrix(unlist(lapply(1:samsize0,function(s){
      res$CDrawsStore[s,g,,][lower.tri(diag(P))]
    })),ncol=samsize0)))
    mean_g <- apply(Fdraws_g,2,mean)
    names(mean_g) <- corrnames[[g]][lower.tri(diag(P))]
    covm_g <- cov(Fdraws_g)
    return(list(mean_g,covm_g))
  })

  meansCovCorr <- lapply(1:numG,function(g){
    matcor_g <- unlist(lapply(1:(P-1),function(p2){
      unlist(lapply((p2+1):P,function(p1){
        mean(res$CDrawsStore[,g,p1,p2])
      }))
    }))
    names(matcor_g) <- corrnames[[g]][lower.tri(diag(P))]
    return(matcor_g)
  })

  meanN <- unlist(lapply(1:numG,function(g){
    FmeansCovCorr[[g]][[1]]
  }))
  covmN <- matrix(0,nrow=numcorr,ncol=numcorr)
  numcorrg <- numcorr/numG

  corrdraws <- lapply(1:numG,function(g){
    array_g <- res$CDrawsStore[,g,,]
    dimnames(array_g) <- list(NULL,varnames[[g]],varnames[[g]])
    return(array_g)
  })

  for(g in 1:numG){
    covmN[(g-1)*numcorrg+1:numcorrg,(g-1)*numcorrg+1:numcorrg] <- FmeansCovCorr[[g]][[2]]
  }

  # posterior estimates
  postestimates_correlations <- Reduce(rbind,
                                       lapply(1:numG,function(g){
                                         means <- meansCovCorr[[g]]
                                         medians <- res$C_quantiles[g,,,2][lower.tri(diag(P))]
                                         lb <- res$C_quantiles[g,,,1][lower.tri(diag(P))]
                                         ub <- res$C_quantiles[g,,,3][lower.tri(diag(P))]
                                         return(cbind(means,medians,lb,ub))
                                       }))
  colnames(postestimates_correlations) <- c("mean","median","2.5%","97.5%")

  cor_out <- list(meanF=meanN,covmF=covmN,correstimates=postestimates_correlations,
                  corrdraws=corrdraws,corrnames=corrnames,variables=varnames)
  class(cor_out) <- "cor_test"

  return(cor_out)
}


