

#' @title Bayesian correlation analysis
#'
#' @name cor_test_old
#'
#' @description Estimate the unconstrained posterior for the correlations using a joint uniform prior.
#'
#' @param ...  matrices (or data frames) of dimensions \emph{n} (observations) by  \emph{p} (variables)
#' for different groups (in case of multiple matrices or data frames).
#'
#' @param formula an object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (e.g., \code{~ education}).
#'
#' @param iter number of iterations from posterior (default is 5000).
#'
#' @return list of class \code{cor_test}:
#' \itemize{
#' \item \code{meanF} posterior means of Fisher transform correlations
#' \item \code{covmF} posterior covariance matrix of Fisher transformed correlations
#' \item \code{correstimates} posterior estimates of correlation coefficients
#' \item \code{corrdraws} list of posterior draws of correlation matrices per group
#' \item \code{corrnames} names of all correlations
#' }
#'
#' @examples
#' \donttest{
#' # Bayesian correlation analysis of the 6 variables in 'memory' object
#' # we consider a correlation analysis of the first three variable of the memory data.
#' fit <- cor_test(BFpack::memory[,1:3])
#'
#' # Bayesian correlation of variables in memory object in BFpack while controlling
#' # for the Cat variable
#' fit <- cor_test(BFpack::memory[,c(1:4)],formula = ~ Cat)
#'
#' # Bayesian correlation analysis of first three variables in memory data
#' # for two different groups
#' HC <- subset(BFpack::memory[,c(1:3,7)], Group == "HC")[,-4]
#' SZ <- subset(BFpack::memory[,c(1:3,7)], Group == "SZ")[,-4]
#' fit <- cor_test(HC,SZ)
#' }
#' @rdname cor_test_old
cor_test_old <- function(..., formula = NULL, iter = 5e3){
  Y_groups <- list(...)

  groups <- length(Y_groups)

  if(is.null(formula)){
    formula <- ~ 1
  }

  model_matrices <- lapply(seq_len(groups) , function(x) {
    model.matrix(formula, Y_groups[[x]])
  })

  correlate <- remove_predictors_helper(Y_groups = Y_groups, formula)

  YXlist <- lapply(1:length(model_matrices),function(g){
    list(as.matrix(correlate[[g]]),as.matrix(model_matrices[[g]]))
  })

  numG <- length(YXlist)
  P <- ncol(YXlist[[1]][[1]])
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
    Y_g <- scale(YXlist[[g]][[1]])
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

  # call Fortran subroutine for Gibbs sampling using noninformative improper priors
  # for regression coefficients, Jeffreys priors for standard deviations, and a proper
  # joint uniform prior for the correlation matrices.
  res <- .Fortran("estimate_postmeancov_fisherz",
                  postZmean=matrix(0,numcorr,1),
                  postZcov=matrix(0,numcorr,numcorr),
                  P=as.integer(P),
                  numcorr=as.integer(numcorr),
                  K=as.integer(K),
                  numG=as.integer(numG),
                  BHat=BHat,
                  sdHat=rbind(sdHat,sdsd),
                  CHat=CHat,
                  XtXi=XtXi,
                  samsize0=as.integer(samsize0),
                  Njs=as.integer(ngroups),
                  Ygroups=Ygroups,
                  Xgroups=Xgroups,
                  Ntot=as.integer(Ntot),
                  C_quantiles=array(0,dim=c(numG,P,P,3)),
                  sigma_quantiles=array(0,dim=c(numG,P,3)),
                  B_quantiles=array(0,dim=c(numG,K,P,3)),
                  BDrawsStore=array(0,dim=c(samsize0,numG,K,P)),
                  sigmaDrawsStore=array(0,dim=c(samsize0,numG,P)),
                  CDrawsStore=array(0,dim=c(samsize0,numG,P,P)),
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









