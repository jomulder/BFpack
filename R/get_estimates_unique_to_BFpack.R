#' @method get_estimates hetcor
#' @import bain
#' @export
get_estimates.hetcor <- function(x, ...){
  cl <- match.call()
  cl[[1]] <- as.name("get_estimates")
  cl[["x"]] <- x$correlations
  P <- nrow(x$std.errors)
  out <- eval.parent(cl)
  retain <- matrix(1:length(out$estimate), nrow = nrow(x$std.errors))
  out$estimate <- out$estimate[retain[lower.tri(retain)]]
  errcov <- x$std.errors**2
  errcov <- errcov[lower.tri(retain)]
  if(length(errcov) == 1){
    out$Sigma <- list(matrix(errcov))
  } else {
    out$Sigma <- list(diag(errcov))
  }
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "hetcor"
  out
}


#' @method get_estimates coxph
#' @export
get_estimates.coxph <- function(x, ...){
  out <- list()
  out$estimate <- coef(x)
  out$Sigma <- list(vcov(x))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "cophx"
  out
}


#' @method get_estimates glm
#' @export
get_estimates.glm <- function(x, ...){
  out <- list()
  out$estimate <- coef(x)
  out$Sigma <- list(vcov(x))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "glm"
  out
}

#' @method get_estimates polr
#' @export
get_estimates.polr <- function(x, ...){
  out <- list()
  out$estimate <- c(coef(x),x$zeta)
  out$Sigma <- list(vcov(x))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "polr"
  out
}


#' @method get_estimates bartlett_htest
#' @export
get_estimates.bartlett_htest <- function(x, ...){
  out <- list()
  out$estimate <- x$vars
  out$Sigma <- NULL
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "bartlett_htest"
  out
}


#' @method get_estimates survreg
#' @export
get_estimates.survreg <- function(x, ...){
  out <- list()
  out$estimate <- x$coefficients
  out$Sigma <- list(x$var)
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "survreg"
  out
}


#' @method get_estimates zeroinfl
#' @export
get_estimates.zeroinfl <- function(x, ...){
  out <- list()
  out$estimate <- c(coef(x),x$zeta)
  out$Sigma <- list(vcov(x))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "zeroinfl"
  out
}


#' @method get_estimates mlm
#' @export
get_estimates.mlm <- function(x, ...){
  out <- list()
  P <- ncol(x$coefficients)
  K <- nrow(x$coefficients)
  N <- nrow(x$residuals)
  names_coef1 <- names(x$coefficients[,1])
  names_coef2 <- names(x$coefficients[1,])
  names_coef <- unlist(lapply(1:P,function(p){
    lapply(1:K,function(k){
      paste0(names_coef1[k],"_on_",names_coef2[p])
    })
  }))
  # estimates of regression coefficients
  estimatesBeta <- c(x$coefficients)
  names(estimatesBeta) <- names_coef
  Xmat <- model.matrix(x)
  Ymat <- model.matrix(x)%*%x$coefficients + x$residuals
  SigmaEst <- t(x$residuals)%*%x$residuals/N
  covmBeta <- kronecker(SigmaEst,solve(t(Xmat)%*%Xmat))
  row.names(covmBeta) <- colnames(covmBeta) <- names_coef

  # estimates of residual correlations
  corrmat <- diag(P)
  numcorrgroup <- P*(P-1)/2
  row.names(corrmat) <- colnames(corrmat) <- colnames(x$residuals)
  corr_names <- names(get_estimates(corrmat)$estimate)
  matrix_names <- matrix(corr_names,nrow=P)
  # equal correlations are at the opposite side of the vector
  corr_names <- c(matrix_names[lower.tri(matrix_names)],
                  t(matrix_names)[lower.tri(matrix_names)])
  CorrMat <- diag(1/sqrt(diag(SigmaEst)))%*%SigmaEst%*%diag(1/sqrt(diag(SigmaEst)))
  estimatesCorr <- rep(CorrMat[lower.tri(CorrMat)],2)
  names(estimatesCorr) <- corr_names

  dummyX <- rep(F,K)
  names(dummyX) <- colnames(Xmat)
  for(k in 1:K){
    uniquek <- sort(unique(Xmat[,k]))
    if(length(uniquek)==2 && uniquek[1]==0 && uniquek[2]==1){dummyX[k]<-T} #group index of intercept
  }
  if(sum(dummyX)==0){
    # no dummy covariates for groups for correlations
    estimatesCorrGroup <- covmCorrGroup <- NULL
    estimates <- c(estimatesBeta,estimatesCorr)
    covm <- diag(rep(1/N,length(estimatesBeta)+length(estimatesCorr)))
    covm[1:length(estimatesBeta),1:length(estimatesBeta)] <- covmBeta
    row.names(covm) <- colnames(covm) <- names(estimates)
  }else{
    # get estimates of correlations per group (identified as dummy covariates)
    corr_estimates_groups <- unlist(lapply(which(dummyX),function(k){
      corr_names_group <- unlist(lapply(1:length(corr_names),function(naam){
        paste0(corr_names[naam],"_in_",names(dummyX[k]))
      }))
      which1 <- which(Xmat[,names(dummyX[k])]==1)
      N1 <- length(which1)
      Xmat1 <- Xmat[which1,]
      Xmat1 <- cbind(rep(1,N1),as.matrix(Xmat1[,apply(Xmat1,2,sd)!=0]))
      Ymat1 <- Ymat[which1,]
      tXX1 <- t(Xmat1)%*%Xmat1
      if(det(tXX1)!=0){
        Bhat1 <- solve(tXX1)%*%t(Xmat1)%*%Ymat1
        SigmaEst1 <- t(Ymat1 - Xmat1%*%Bhat1)%*%(Ymat1 - Xmat1%*%Bhat1)/N1
        CorrMat1 <- diag(1/sqrt(diag(SigmaEst1)))%*%SigmaEst1%*%
          diag(1/sqrt(diag(SigmaEst1)))
        estimatesCorr1 <- rep(CorrMat[lower.tri(CorrMat1)],2)
        names(estimatesCorr1) <- corr_names_group
        return(list(estimatesCorr1))
      }
    }))
    Ngroups <- apply(Xmat[,dummyX],2,sum)
    #combine estimates and covariance matrix
    estimates <- c(estimatesBeta,estimatesCorr,corr_estimates_groups)
    covm <- diag(c(rep(0,length(estimatesBeta)),rep(1/N,length(estimatesCorr)),rep(1/Ngroups,each=length(estimatesCorr))))
    covm[1:length(estimatesBeta),1:length(estimatesBeta)] <- covmBeta
    row.names(covm) <- colnames(covm) <- names(estimates)
  }
  out$estimate <- estimates
  out$Sigma <- covm
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "mlm"
  out
}

