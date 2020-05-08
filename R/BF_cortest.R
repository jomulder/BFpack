
#' @title Bayesian correlation analysis
#'
#' @name cor_test
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
#'
#' # Bayesian correlation analysis of the 6 variables in 'memory' object
#' # sample posterior
#' fit <- cor_test(BFpack::memory[,1:6])
#'
#' # Bayesian correlation of variables in memory object in BFpack while controlling
#' # for the Group variable
#' # sample posterior
#' fit <- cor_test(BFpack::memory,formula = ~ Group)
#'
#' # Bayesian correlation of variables in memory object in BFpack separate for the
#' # two different groups
#' # sample posterior
#' HC <- subset(BFpack::memory, Group == "HC")[,-7]
#' SZ <- subset(BFpack::memory, Group == "SZ")[,-7]
#' fit <- cor_test(HC,SZ)
#'
#' @rdname cor_test
#' @export
cor_test <- function(..., formula = NULL, iter = 5e3){
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
  }
  samsize0 <- iter

  # call Fortran subroutine for Gibbs sampling using noninformative improper priors
  # for regression coefficients, Jeffreys priors for standard deviations, and a proper
  # joint uniform prior for the correlation matrices.
  res <-.Fortran("estimate_postmeancov_fisherz",
                 postZmean=matrix(0,numcorr,1),
                 postZcov=matrix(0,numcorr,numcorr),
                 P=as.integer(P),
                 numcorr=as.integer(numcorr),
                 K=as.integer(K),
                 numG=as.integer(numG),
                 BHat=BHat,
                 sdHat=sdHat,
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
                  corrdraws=corrdraws,corrnames=corrnames)
  class(cor_out) <- "cor_test"

  return(cor_out)
}

#' @importFrom stats terms
remove_predictors_helper <- function(Y_groups, formula){

  # number of groups
  groups <- length(Y_groups)

  # model matrix terms
  mm_terms <- attr(terms(formula), "term.labels")

  if(length(mm_terms) == 0){

    Y_groups

  } else {
    lapply(seq_len(groups), function(x){

      # check for factors
      factor_pred <- which(paste0("as.factor(", colnames(Y_groups[[x]]), ")") %in% mm_terms)

      # check for non factors
      cont_pred <- which(colnames(Y_groups[[x]]) %in% mm_terms)

      # remove predictors
      Y_groups[[x]][,-c(factor_pred, cont_pred)]

    })
  }
}

FisherZ <- function(r){.5*log((1+r)/(1-r))}


#' @importFrom mvtnorm dmvnorm pmvnorm rmvnorm
#' @importFrom stats dnorm pnorm
#' @method BF cor_test
#' @export
BF.cor_test <- function(x,
                        hypothesis = NULL,
                        prior = NULL,
                        ...){

  bayesfactor <- "Bayes factor based on joint uniform priors for correlations"
  testedparameter <- "correlation coefficients"

  P <- dim(x$corrdraws[[1]])[2]
  numG <- length(x$corrdraws)
  numcorrgroup <- P*(P-1)/2
  get_est <- get_estimates(x)
  corrmeanN <- get_est$estimate
  corrcovmN <- get_est$Sigma[[1]]

  # Exploratory testing of correlation coefficients
  #get height of prior density at 0 of Fisher transformed correlation
  drawsJU <- draw_ju_r(P,samsize=50000,Fisher=1)
  relcomp0 <- approxfun(density(drawsJU[,1]))(0)
  # compute exploratory BFs
  corr_names <- rownames(x$correstimates)
  numcorr <- length(corrmeanN)
  relfit <- matrix(c(dnorm(0,mean=corrmeanN,sd=sqrt(diag(corrcovmN))),
                     pnorm(0,mean=corrmeanN,sd=sqrt(diag(corrcovmN))),
                     1-pnorm(0,mean=corrmeanN,sd=sqrt(diag(corrcovmN)))),ncol=3)
  relcomp <- matrix(c(rep(relcomp0,numcorr),rep(.5,numcorr*2)),ncol=3)
  colnames(relcomp) <- colnames(relfit) <- c("p(=0)","Pr(<0)","Pr(>0)")
  BFtu_exploratory <- relfit / relcomp
  row.names(BFtu_exploratory) <- rownames(x$correstimates)
  colnames(BFtu_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
  PHP_exploratory <- round(BFtu_exploratory /
                             apply(BFtu_exploratory,1,sum),3)
  # posterior estimates
  postestimates <- x$correstimates

  # confirmatory testing if hypothesis argument is used
  if(!is.null(hypothesis)){

    #check if constraints are formulated on correlaties in different populations
    #if so, then the correlation names contains the string "_group" at the end
    params_in_hyp1 <- params_in_hyp(hypothesis)

    corr_names <- unlist(lapply(1:length(x$corrnames),function(g){
      c(x$corrnames[[g]][lower.tri(x$corrnames[[g]])],
        t(x$corrnames[[g]])[lower.tri(x$corrnames[[g]])])
    })) #which includes Y1_with_Y2 and Y2_with_Y1
    # checkLabels <- matrix(unlist(lapply(1:length(params_in_hyp1),function(par){
    #   params_in_hyp1[par]==corr_names
    # })),nrow=length(params_in_hyp1),byrow=T)

    parse_hyp <- parse_hypothesis(corr_names,hypothesis)
    parse_hyp$hyp_mat <- do.call(rbind, parse_hyp$hyp_mat)
    if(nrow(parse_hyp$hyp_mat)==1){
      select1 <- rep(1:numcorrgroup,numG) + rep((0:(numG-1))*2*numcorrgroup,each=numcorrgroup)
      select2 <- rep(numcorrgroup+1:numcorrgroup,numG) + rep((0:(numG-1))*2*numcorrgroup,each=numcorrgroup)
      parse_hyp$hyp_mat <-
        t(as.matrix(c(parse_hyp$hyp_mat[,select1] + parse_hyp$hyp_mat[,select2],parse_hyp$hyp_mat[,numcorrgroup*2*numG+1])))
    }else{
      #combine equivalent correlations, e.g., cor(Y1,Y2)=corr(Y2,Y1).
      select1 <- rep(1:numcorrgroup,numG) + rep((0:(numG-1))*2*numcorrgroup,each=numcorrgroup)
      select2 <- rep(numcorrgroup+1:numcorrgroup,numG) + rep((0:(numG-1))*2*numcorrgroup,each=numcorrgroup)
      parse_hyp$hyp_mat <-
        cbind(parse_hyp$hyp_mat[,select1] + parse_hyp$hyp_mat[,select2],parse_hyp$hyp_mat[,numcorrgroup*2*numG+1])
    }
    #create coefficient with equality and order constraints
    RrList <- make_RrList2(parse_hyp)
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]

    # RrStack <- rbind(do.call(rbind,RrE),do.call(rbind,RrO))
    # RStack <- RrStack[,-(numcorr+1)]
    # rStack <- RrStack[,(numcorr+1)]

    numhyp <- length(RrE)
    relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Gaussian_measures(corrmeanN,corrcovmN,RrE1=RrE[[h]],RrO1=RrO[[h]])
    })),nrow=2))
    relcomp <- t(matrix(unlist(lapply(1:numhyp,function(h){
      jointuniform_measures(P,numcorrgroup,numG,RrE1=RrE[[h]],RrO1=RrO[[h]],Fisher=1)
    })),nrow=2))
    row.names(relfit) <- row.names(relcomp) <- parse_hyp$original_hypothesis
    relfit <- Gaussian_prob_Hc(corrmeanN,corrcovmN,relfit,RrO)
    relcomp <- jointuniform_prob_Hc(P,numcorrgroup,numG,relcomp,RrO)
    hypothesisshort <- unlist(lapply(1:nrow(relfit),function(h) paste0("H",as.character(h))))
    row.names(relfit) <- row.names(relfit) <- hypothesisshort

    # the BF for the complement hypothesis vs Hu needs to be computed.
    BFtu_confirmatory <- c(apply(relfit / relcomp, 1, prod))
    # Check input of prior probabilies
    if(is.null(prior)){
      priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
    }else{
      if(!is.numeric(prior) || length(prior)!=length(BFtu_confirmatory)){
        warning(paste0("Argument 'prior' should be numeric and of length ",as.character(length(BFtu_confirmatory)),". Equal prior probabilities are used."))
        priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
      }else{
        priorprobs <- prior
      }
    }
    names(priorprobs) <- names(BFtu_confirmatory)
    PHP_confirmatory <- BFtu_confirmatory*priorprobs / sum(BFtu_confirmatory*priorprobs)
    BFtable <- cbind(relcomp,relfit,relfit[,1]/relcomp[,1],relfit[,2]/relcomp[,2],
                     apply(relfit,1,prod)/apply(relcomp,1,prod),PHP_confirmatory)
    row.names(BFtable) <- names(BFtu_confirmatory)
    colnames(BFtable) <- c("complex=","complex>","fit=","fit>","BF=","BF>","BF","PHP")
    BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory))/
      t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
    row.names(BFmatrix_confirmatory) <- colnames(BFmatrix_confirmatory) <- names(BFtu_confirmatory)
    if(nrow(relfit)==length(parse_hyp$original_hypothesis)){
      hypotheses <- parse_hyp$original_hypothesis
    }else{
      hypotheses <- c(parse_hyp$original_hypothesis,"complement")
    }

  }else{
    BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relfit <-
      relcomp <- hypotheses <- BFtable <- priorprobs <- NULL
  }

  BFcorr_out <- list(
    BFtu_exploratory=BFtu_exploratory,
    PHP_exploratory=PHP_exploratory,
    BFtu_confirmatory=BFtu_confirmatory,
    PHP_confirmatory=PHP_confirmatory,
    BFmatrix_confirmatory=BFmatrix_confirmatory,
    BFtable_confirmatory=BFtable,
    prior=priorprobs,
    hypotheses=hypotheses,
    estimates=postestimates,
    model=x,
    bayesfactor=bayesfactor,
    parameter=testedparameter,
    call=match.call())

  class(BFcorr_out) <- "BF"

  return(BFcorr_out)

}


