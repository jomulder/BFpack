#' @importFrom pracma rref
#' @importFrom mvtnorm dmvnorm pmvnorm rmvnorm dmvt pmvt
#' @importFrom Matrix rankMatrix
#' @importFrom MCMCpack rinvgamma
#' @importFrom MASS ginv
#' @method BF mlm
#' @export
BF.mlm <- function(x,
                   hypothesis = NULL,
                   prior = NULL,
                   ...){

  # default BF on location parameters in a univarite normal linear model
  # Note that it is recommended that the fitten model is based on standardized covariates.
  P <- ncol(x$residuals)
  N <- nrow(x$residuals)
  K <- length(x$coefficients)/P # dimension of predictors per dependent variable
  dummyX <- rep(F,K)
  names(dummyX) <- row.names(x$coefficients)

  Xmat <- model.matrix(x)
  Ymat <- model.matrix(x)%*%x$coefficients + x$residuals

  # Exploratory testing of regression coefficients
  if(length(x$xlevels)==0){ #no grouping covariates: 1 group
    J <- 1
    dummyX <- rep(F,K)
    Nj <- nrow(Xmat)
    dvec <- rep(1,Nj)
    #set minimal fractions for the group
    bj <- ((P+K)/J)/Nj
  }else{
    numlevels <- unlist(lapply(x$xlevels,length))
    mains <- unlist(lapply(1:length(x$xlevels),function(fac){
      unlist(lapply(1:length(x$xlevels[[fac]]),function(lev){
        paste0(names(x$xlevels)[fac],x$xlevels[[fac]][lev])
      }))
    }))
    intercept <- attr(x$terms,"intercept")==1
    names_coef <- row.names(x$coefficients)
    # dummyX indicate which columns contain dummy group covariates
    dummyX1 <- apply(matrix(unlist(lapply(1:length(mains),function(faclev){
      unlist(lapply(1:length(names_coef),function(cf){
        grepl(mains[faclev],names_coef[cf])
      }))
    })),nrow=length(names_coef)),1,max)==1

    dummyX2 <- unlist(lapply(1:ncol(Xmat),function(coldum){
      length(table(Xmat[,coldum]))
    }))==2
    dummyX <- dummyX1 * dummyX2 == 1
    #number of groups on variations of dummy combinations
    groupcode <- as.matrix(unique(Xmat[,dummyX]))
    rownames(groupcode) <- unlist(lapply(1:nrow(groupcode),function(r){
      paste0("groupcode",r)
    }))
    J <- nrow(groupcode)
    if(J==nrow(Xmat)){
      stop("Not enough observations for every group. Try fitting the model without factors.")
    }

    # group membership of each observation
    dvec <- unlist(lapply(1:N,function(i){
      which(rowSums(abs(t(matrix(rep(Xmat[i,dummyX],J),ncol=J)) - groupcode))==0)
    }))
    Nj <- c(table(dvec))
    #set minimal fractions for each group
    bj <- ((P+K)/J)/Nj
  }

  #Compute sufficient statistics for all groups
  tXXj <- lapply(1:J,function(j){
    if(Nj[j]==1){
      Xmat[dvec==j,]%*%t(Xmat[dvec==j,])
    }else t(Xmat[dvec==j,])%*%Xmat[dvec==j,]
  })
  tXXj_b <- lapply(1:J,function(j){
    tXXj[[j]]*bj[j]
  })
  tXYj <- lapply(1:J,function(j){
    if(Nj[j]==1){
      as.matrix(Xmat[dvec==j,]*Ymat[dvec==j,])
    } else {t(Xmat[dvec==j,])%*%Ymat[dvec==j,]}
  })
  tXYj_b <- lapply(1:J,function(j){
    tXYj[[j]]*bj[j]
  })
  tYYj <- lapply(1:J,function(j){
    t(Ymat[dvec==j,])%*%Ymat[dvec==j,]
  })
  tYYj_b <- lapply(1:J,function(j){
    tYYj[[j]]*bj[j]
  })
  tXX <- Reduce("+",tXXj)
  tXXi <- solve(tXX)
  tXY <- Reduce("+",tXYj)
  tYY <- Reduce("+",tYYj)
  tXX_b <- Reduce("+",tXXj_b)
  tXXi_b <- solve(tXX_b)
  tXY_b <- Reduce("+",tXYj_b)
  tYY_b <- Reduce("+",tYYj_b)
  BetaHat <- solve(tXX)%*%tXY          # same as x$coefficients
  S <- tYY - t(tXY)%*%solve(tXX)%*%tXY # same as sum((x$residuals)**2)
  # sufficient statistics based on fraction of the data
  BetaHat_b <- solve(tXX_b)%*%tXY_b
  S_b <- tYY_b - t(tXY_b)%*%solve(tXX_b)%*%tXY_b

  # BF computation for exploratory analysis of separate parameters
  names_coef1 <- names(x$coefficients[,1])
  names_coef2 <- names(x$coefficients[1,])
  names_coef <- unlist(lapply(1:P,function(p){
    lapply(1:K,function(k){
      paste0(names_coef1[k],"_on_",names_coef2[p])
    })
  }))

  # prior hyperparameters
  df0 <- 1 # should be the same as sum(rep(bj,times=Nj))-K-P+1
  Scale0 <- kronecker(S_b,tXXi_b)
  mean0 <- as.matrix(rep(0,K*P))
  # posterior hyperparameters
  dfN <- N-K-P+1
  ScaleN <- kronecker(S,tXXi)/(N-K-P+1) # off-diagonal elements have no meaning
  meanN <- as.matrix(c(BetaHat))
  row.names(meanN) <- row.names(mean0) <- names_coef

  # Hypotheses for exploratory test
  # H0: beta = 0
  # H1: beta < 0
  # H2: beta < 0
  relfit <- t(matrix(unlist(lapply(1:(K*P),function(k){
    c(dt((0-meanN[k,1])/sqrt(ScaleN[k,k]),df=dfN)/sqrt(ScaleN[k,k]),
      pt((0-meanN[k,1])/sqrt(ScaleN[k,k]),df=dfN,lower.tail = TRUE),
      pt((0-meanN[k,1])/sqrt(ScaleN[k,k]),df=dfN,lower.tail = FALSE))
  })),nrow=3))
  relcomp <- t(matrix(unlist(lapply(1:(K*P),function(k){
    c(dt((0-mean0[k,1])/sqrt(Scale0[k,k]),df=df0)/sqrt(Scale0[k,k]),
      pt((0-mean0[k,1])/sqrt(Scale0[k,k]),df=df0,lower.tail = TRUE),
      pt((0-mean0[k,1])/sqrt(Scale0[k,k]),df=df0,lower.tail = FALSE))
  })),nrow=3))
  colnames(relfit) <- colnames(relcomp) <- c("p(=0)","Pr(<0)","Pr(>0)")
  row.names(relcomp) <- row.names(relfit) <- names_coef

  BFtu_exploratory_coefficients <- relfit / relcomp
  colnames(BFtu_exploratory_coefficients) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
  PHP_exploratory_coefficients <- BFtu_exploratory_coefficients /
    apply(BFtu_exploratory_coefficients,1,sum)

  #compute estimates
  postestimates_coefficients <- cbind(meanN,meanN,
                                      t(matrix(unlist(lapply(1:length(meanN),function(coef){
                                        ub <- qt(p=.975,df=dfN)*sqrt(ScaleN[coef,coef])+meanN[coef,1]
                                        lb <- qt(p=.025,df=dfN)*sqrt(ScaleN[coef,coef])+meanN[coef,1]
                                        return(c(ub,lb))
                                      })),nrow=2))
  )
  row.names(postestimates_coefficients) <- names_coef
  colnames(postestimates_coefficients) <- c("mean","median","2.5%","97.5%")

  # testedparameter <- "regression coefficients"

  # Additional exploratory tests of main effects and interaction effects
  # in the case of an aov type object
  if(sum(class(x)=="aov")==1){
    testedparameter <- "group means"

    # check main effects
    BFmain <- unlist(lapply(1:length(numlevels),function(fac){
      name1 <- names(numlevels[fac])
      mains1 <- mains[sum(numlevels[1:fac])-numlevels[fac]+1:numlevels[fac]]
      which0 <- unlist(lapply(1:length(colnames(Xmat)),function(col){
        sum(colnames(Xmat)[col]==mains1)==1
      }))
      if(sum(which0)>0){
        RrE_f <- matrix(0,nrow=sum(which0),ncol=length(colnames(Xmat)))
        for(r1 in 1:sum(which0)){RrE_f[r1,which(which0)[r1]]<-1}
        RrE_f <- cbind(kronecker(diag(P),RrE_f),rep(0,sum(which0)*P))
        relcomp_f <- MatrixStudent_measures(Mean1=matrix(mean0,ncol=P),Scale1=S_b,tXXi1=tXXi_b,
                                            df1=df0,RrE1=RrE_f,RrO1=NULL,Names1=NULL,constraints1=NULL,
                                            MCdraws=1e4)
        relfit_f <- MatrixStudent_measures(Mean1=BetaHat,Scale1=S,tXXi1=tXXi,df1=dfN,RrE1=RrE_f,
                                           RrO1=NULL,Names1=NULL,constraints1=NULL,MCdraws=1e4)
        BFtu <- relfit_f[1]/relcomp_f[1]
        names(BFtu) <- name1
        return(c(BFtu,relfit_f[1],relcomp_f[1]))
      }
    }))
    #compute Bayes factors for testing main effects if present
    if(length(BFmain)>0){ # then there are main effects
      names_main <- names(BFmain[(0:(length(BFmain)/3-1))*3+1])
      BFtu_main <- matrix(c(BFmain[(0:(length(BFmain)/3-1))*3+1],rep(1,length(BFmain)/3)),nrow=length(BFmain)/3)
      row.names(BFtu_main) <- names_main
      colnames(BFtu_main) <- c("BFtu","BFuu")
      PHP_main <- BFtu_main / apply(BFtu_main,1,sum)
      colnames(PHP_main) <- c("Pr(null)","Pr(alt)")
    }else{ PHP_main <- BFtu_main <- NULL}
    #check whether interaction effects are present
    prednames <- names(attr(x$term,"dataClasses"))
    matcov <- cbind(matrix(unlist(lapply(1:K,function(col){
      colx <- colnames(Xmat)[col]
      unlist(lapply(1:length(prednames),function(pred){
        grepl(prednames[pred],colx)
      }))
    })),nrow=length(prednames)),rep(F,length(prednames)))
    row.names(matcov) <- prednames
    colnames(matcov) <- c(colnames(Xmat),"dummy")
    BFtu_interaction0 <- list()
    count_interaction <- 0

    # cat("Computing Bayes factors for testing interaction effects may take a few seconds.")
    # cat("\n")
    # cat("\n")

    for(c1 in 1:ncol(matcov)){
      if(c1 < ncol(matcov)){
        numeffects_c <- sum(matcov[,c1])
        if(numeffects_c>1){
          count_interaction <- count_interaction + 1
          interactionset <- names(which(matcov[,c1]))
          whichx <- apply(matcov[which(matcov[,c1]),],2,sum)==length(interactionset)
          RrE_ia <- matrix(0,nrow=sum(whichx),ncol=K)
          for(r1 in 1:sum(whichx)){RrE_ia[r1,which(whichx)[r1]]<-1}
          RrE_ia <- cbind(kronecker(diag(P),RrE_ia),rep(0,sum(whichx)*P))
          relcomp_ia <- MatrixStudent_measures(Mean1=matrix(mean0,ncol=P),Scale1=S_b,tXXi1=tXXi_b,
                                               df1=df0,RrE1=RrE_ia,RrO1=NULL,Names1=NULL,
                                               constraints1=NULL,MCdraws=1e4)
          names(relcomp_ia) <- c("c=","c>")
          relfit_ia <- MatrixStudent_measures(Mean1=BetaHat,Scale1=S,tXXi1=tXXi,df1=dfN,RrE1=RrE_ia,
                                              RrO1=NULL,Names1=NULL,constraints1=NULL,MCdraws=1e4)
          names(relfit_ia) <- c("f=","f>")
          BFtu_ia <- relfit_ia[1]/relcomp_ia[1]
          names(BFtu_ia) <- paste(interactionset,collapse=":")
          BFtu_interaction0[[count_interaction]] <- c(BFtu_ia,relcomp_ia[1],relfit_ia[1])
          #exclude other columns from Xmat that have been used to avoid double results
          matcov <- matcov[,-(which(whichx[-c1])+1)]
        }
      }
    }
    #compute Bayes factors for testing interaction effects if present
    if(count_interaction>0){ # then there are main effects
      BFtu_interaction0 <- unlist(BFtu_interaction0)
      names_interaction <- names(BFtu_interaction0[(0:(length(BFtu_interaction0)/3-1))*3+1])
      BFtu_interaction <- matrix(c(BFtu_interaction0[(0:(length(BFtu_interaction0)/3-1))*3+1],
                                   rep(1,length(BFtu_interaction0)/3)),nrow=length(BFtu_interaction0)/3)
      row.names(BFtu_interaction) <- names_interaction
      colnames(BFtu_interaction) <- c("BFtu","BFuu")
      PHP_interaction <- BFtu_interaction / apply(BFtu_interaction,1,sum)
      colnames(PHP_interaction) <- c("Pr(null)","Pr(alt)")
    }else{ PHP_interaction <- BFtu_interaction <- NULL}
  }else{ PHP_interaction <- BFtu_interaction <- PHP_main <- BFtu_main <- NULL}

  # Exploratory testing of correlation coefficients
  corrmat <- diag(P)
  numcorrgroup <- P*(P-1)/2
  row.names(corrmat) <- colnames(corrmat) <- colnames(x$residuals)
  corr_names <- names(get_estimates(corrmat)$estimate)
  matrix_names <- matrix(corr_names,nrow=P)
  # equal correlations are at the opposite side of the vector
  corr_names <- c(matrix_names[lower.tri(matrix_names)],
                  t(matrix_names)[lower.tri(matrix_names)])
  #exploratory tests on the residual correlations of the mlm object
  numG <- 1
  ngroups <- nrow(Ymat)
  YXlist <- list(list(Ymat,Xmat))
  corr_names_exploratory <- matrix_names[lower.tri(matrix_names)]

  # Gibbs sampler to estimate the posterior of the Fisher transformed
  # correlations having an approximate Gaussian distribution.
  Gibbs_output <- estimate_postMeanCov_FisherZ(YXlist,numdraws=8e3)
  corrmeanN <- Gibbs_output$meanN
  corrcovmN <- Gibbs_output$covmN
  numcorr <- length(corrmeanN)

  #get height of prior density at 0 of Fisher transformed correlation
  drawsJU <- draw_ju_r(P,samsize=50000,Fisher=1)
  relcomp0 <- approxfun(density(drawsJU[,1]))(0)
  # compute exploratory BFs
  relfit <- matrix(c(dnorm(0,mean=corrmeanN,sd=sqrt(diag(corrcovmN))),
                     pnorm(0,mean=corrmeanN,sd=sqrt(diag(corrcovmN))),
                     1-pnorm(0,mean=corrmeanN,sd=sqrt(diag(corrcovmN)))),ncol=3)
  relcomp <- matrix(c(rep(relcomp0,numcorr),rep(.5,numcorr*2)),ncol=3)
  colnames(relcomp) <- colnames(relfit) <- c("p(=0)","Pr(<0)","Pr(>0)")
  BFtu_exploratory_correlations <- relfit / relcomp
  row.names(BFtu_exploratory_correlations) <- corr_names_exploratory
  colnames(BFtu_exploratory_correlations) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
  PHP_exploratory_correlations <- round(BFtu_exploratory_correlations /
                                          apply(BFtu_exploratory_correlations,1,sum),3)
  # posterior estimates
  postestimates_correlations <- Reduce(rbind,
                                       lapply(1:numG,function(g){
                                         means <- Gibbs_output$corr_means[[g]]
                                         medians <- Gibbs_output$corr_quantiles[g,,,2][lower.tri(diag(P))]
                                         lb <- Gibbs_output$corr_quantiles[g,,,1][lower.tri(diag(P))]
                                         ub <- Gibbs_output$corr_quantiles[g,,,3][lower.tri(diag(P))]
                                         return(cbind(means,medians,lb,ub))
                                       }))
  row.names(postestimates_correlations) <- corr_names_exploratory
  colnames(postestimates_correlations) <- c("mean","median","2.5%","97.5%")

  #combine BFs and PHPs of exploratory tests of coefficients and correlations
  BFtu_exploratory <- rbind(BFtu_exploratory_coefficients,BFtu_exploratory_correlations)
  PHP_exploratory <- rbind(PHP_exploratory_coefficients,PHP_exploratory_correlations)
  postestimates <- rbind(postestimates_coefficients,postestimates_correlations)

  # confirmatory BF test
  if(!is.null(hypothesis)){

    # hypotheses are either formulated only on correlations, or only on regression
    # coefficients
    if(grepl("_with_",hypothesis)){
      #then constraints on correlations

      #check if constraints are formulated on correlaties in different populations
      #if so, then the correlation names contains the string "_group" at the end
      params_in_hyp1 <- params_in_hyp(hypothesis)
      # check if labels of correlations in hypotheses match with actual correlations
      checkLabels <- matrix(unlist(lapply(1:length(params_in_hyp1),function(par){
        params_in_hyp1[par]==corr_names
      })),nrow=length(params_in_hyp1),byrow=T)
      if(sum(abs(apply(checkLabels,1,sum)-rep(1,length(params_in_hyp1))))==0){
        #constraints are formulated on correlations from one population/group
        numG <- 1
        ngroups <- nrow(Ymat)
        YXlist <- list(list(Ymat,Xmat))
      }else{
        #check if labels of correlations in hypotheses match with
        #combined labels of correlations and dummy group variables

        #Check which are dummy variables corresponding to (adjusted) mean parameters
        dummyX <- rep(F,K)
        names(dummyX) <- colnames(Xmat)
        for(k in 1:K){
          uniquek <- sort(unique(Xmat[,k]))
          if(length(uniquek)==2 && uniquek[1]==0 && uniquek[2]==1){dummyX[k]<-T} #group index of intercept
        }
        if(sum(dummyX)==0){
          # there appears to be no dummy group variables
          stop("Labels for correlations should be of the form 'y1_with_y2' or 'y1_with_y2_in_x1', where x1 must be a dummy (0/1) group covariate.")
        }
        # create correlation labels
        corr_names_all <- unlist(lapply(which(dummyX==T),function(k){
          if(dummyX[k]){
            unlist(lapply(1:length(corr_names),function(naam){
              paste0(corr_names[naam],"_in_",names(dummyX[k]))
            }))
          }
        }))
        names(corr_names_all) <- NULL
        checkLabels <- matrix(unlist(lapply(1:length(params_in_hyp1),function(par){
          params_in_hyp1[par]==corr_names_all
        })),nrow=length(params_in_hyp1),byrow=T)

        if(sum(abs(apply(checkLabels,1,sum)-rep(1,length(params_in_hyp1))))!=0){
          stop("Labels for correlations should be of the form 'y1_with_y2' or 'y1_with_y2_in_x1', where x1 must be a dummy (0/1) group covariate.")
        }
        #find which covariates are group specific for the correlations on which hypotheses are formulated.
        groupcorrelations <- corr_names_all[which(apply(checkLabels,2,sum)==1)]
        groupcov <- unique(unlist(lapply(1:length(groupcorrelations),function(grcor){
          substr(groupcorrelations[grcor],start=max(gregexpr(pattern ='_',groupcorrelations[grcor])[[1]])+1,stop=nchar(groupcorrelations[grcor]))
        })))
        Xgroup <- matrix(unlist(lapply(1:length(groupcov),function(gcov){
          dvec <- Xmat[,groupcov[gcov]]
        })),nrow=nrow(Xmat))
        colnum <- unlist(lapply(1:length(groupcov),function(gcov){
          which(colnames(Xmat)==groupcov[gcov])
        }))
        colnames(Xgroup) <- groupcov
        which0 <- which(apply(Xgroup,1,sum)==0)
        if(length(which0)>0){
          grouprest <- rep(0,nrow(Xgroup))
          grouprest[which0] <- 1
          Xgroup <- cbind(Xgroup,grouprest)
          colnames(Xgroup) <- c(groupcov,"restgroup")
        }
        ngroups <- apply(Xgroup,2,sum)
        numG <- ncol(Xgroup) # numG
        whichgr1 <- which(apply(Xgroup,1,sum)>1)
        if(length(whichgr1)>0){ #extra check for dummy group variables
          stop("Group indices for dummy group variables need to be either 0 or 1.")
        }
        YXlist <- lapply(1:numG,function(col){
          return(list(Y=Ymat[Xgroup[,col]==1,],
                      X=as.matrix(Xmat[Xgroup[,col]==1,-colnum[-col]])))
        })
        corr_names <- unlist(lapply(1:length(groupcov),function(gr){
          unlist(lapply(1:length(corr_names),function(cn){
            paste0(corr_names[cn],"_in_",groupcov[gr])
          }))
        }))
      }
      # if correlations are tested in different groups get posterior of
      # Fisher transformed correlations per group
      if(numG>1){
        Gibbs_output <- estimate_postMeanCov_FisherZ(YXlist,numdraws=8e3)
        corrmeanN <- Gibbs_output$meanN
        corrcovmN <- Gibbs_output$covmN
        numcorr <- length(corrmeanN)
      }

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

      RrStack <- rbind(do.call(rbind,RrE),do.call(rbind,RrO))
      RStack <- RrStack[,-(numcorr+1)]
      rStack <- RrStack[,(numcorr+1)]

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
      colnames(BFtable) <- c("comp_E","comp_O","fit_E","fit_O","BF_E","BF_O","BF","PHP")
      BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory))/
        t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
      row.names(BFmatrix_confirmatory) <- colnames(BFmatrix_confirmatory) <- names(BFtu_confirmatory)
      if(nrow(relfit)==length(parse_hyp$original_hypothesis)){
        hypotheses <- parse_hyp$original_hypothesis
      }else{
        hypotheses <- c(parse_hyp$original_hypothesis,"complement")
      }

    }else{
      #then constraints on regression coefficients

      #read constraints
      names_coef1 <- names(x$coefficients[,1])
      names_coef2 <- names(x$coefficients[1,])
      names_coef <- unlist(lapply(1:P,function(p){
        lapply(1:K,function(k){
          paste(names_coef1[k],"_on_",names_coef2[p],sep="")
        })
      }))
      matrixnames <- matrix(names_coef,nrow=K)

      # translate named constraints to matrices with coefficients for constraints
      parse_hyp <- parse_hypothesis(names_coef,hypothesis) #, return_list = TRUE)
      parse_hyp$hyp_mat <- do.call(rbind, parse_hyp$hyp_mat)
      RrList <- make_RrList2(parse_hyp)
      RrE <- RrList[[1]]
      RrO <- RrList[[2]]
      # row.names(RrO[[1]]) <- colnames(RrO[[1]]) <- row.names(RrE[[1]]) <- colnames(RrE[[1]]) <- NULL
      # row.names(RrO[[2]]) <- colnames(RrO[[2]]) <- row.names(RrE[[2]]) <- colnames(RrE[[2]]) <- NULL

      RrStack <- rbind(do.call(rbind,RrE),do.call(rbind,RrO))
      RrStack <- interval_RrStack(RrStack)
      if(nrow(RrStack)>1){
        RStack <- RrStack[,-(K*P+1)]
        rStack <- RrStack[,(K*P+1)]
      }else{
        RStack <- matrix(RrStack[,-(K*P+1)],nrow=1)
        rStack <- RrStack[,(K*P+1)]
      }

      # check if a common boundary exists for prior location under all constrained hypotheses
      if(nrow(RrStack) > 1){
        rref_ei <- rref(RrStack)
        nonzero <- RrStack[,K+1]!=0
        if(max(nonzero)>0){
          row1 <- max(which(nonzero==T))
          if(sum(abs(rref_ei[row1,1:K]))==0){
            stop("No common boundary point for prior location. Conflicting constraints.")
          }
        }
      }

      #number of hypotheses that are specified
      numhyp <- length(RrO)
      #default prior location
      Mean0 <- matrix(c(ginv(RStack)%*%rStack),nrow=K,ncol=P)

      relmeasunlist <- unlist(lapply(1:numhyp,function(h){
        # Check whether the constraints are on a single row or column, if so
        # use the analytic expression, else using a Monte Carlo estimate.
        RrStack <- rbind(RrE[[h]],RrO[[h]])
        Rcheck <- Reduce("+",lapply(1:nrow(RrStack),function(row1){
          abs(matrix(RrStack[row1,-(K*P+1)],nrow=K))
        }))
        RcheckRow <- apply(Rcheck,1,sum)
        RcheckCol <- apply(Rcheck,2,sum)
        if(sum(RcheckRow!=0)==1){ # use multivariate Student distributions
          K1 <- which(RcheckRow!=0)
          # posterior hyperparameters
          dfN <- N-K-P+1
          ScaleN <- S*tXXi[K1,K1]/(N-K-P+1) # off-diagonal elements have no meaning
          meanN <- as.matrix(c(BetaHat[K1,]))
          # exclude inactive rows
          if(is.null(RrE[[h]])){RrE_h=NULL
          }else{
            if(nrow(RrE[[h]])==1){
              RrE_h <- t(as.matrix(RrE[[h]][,c((0:(P-1))*K+K1,P*K+1)]))
            }else{
              RrE_h <- RrE[[h]][,c((0:(P-1))*K+K1,P*K+1)]
            }
          }
          if(is.null(RrO[[h]])){RrO_h=NULL
          }else{
            if(nrow(RrO[[h]])==1){
              RrO_h <- t(as.matrix(RrO[[h]][,c((0:(P-1))*K+K1,P*K+1)]))
            }else{
              RrO_h <- RrO[[h]][,c((0:(P-1))*K+K1,P*K+1)]
            }
          }
          # prior hyperparameters
          df0 <- 1 # should be the same as sum(rep(bj,times=Nj))-K-P+1
          Scale0 <- S_b*tXXi_b[K1,K1]
          mean0 <- Mean0[K1,]
          # compute relative measures of fit and complexity
          relcomp_h <- Student_measures(mean1=mean0,Scale1=Scale0,df1=df0,RrE1=RrE_h,
                                        RrO1=RrO_h,names1=matrixnames[K1,],
                                        constraints1=parse_hyp$original_hypothesis[h])
          relfit_h <- Student_measures(mean1=meanN,Scale1=ScaleN,df1=dfN,RrE1=RrE_h,
                                       RrO1=RrO_h,names1=matrixnames[K1,],
                                       constraints1=parse_hyp$original_hypothesis[h])

        }else if(sum(RcheckCol!=0)==1){ # use multivariate Student distributions
          P1 <- which(RcheckCol!=0)
          # prior hyperparameters
          df0 <- 1 # should be the same as sum(rep(bj,times=Nj))-K-P+1
          Scale0 <- S_b[P1,P1]*tXXi_b
          mean0 <- Mean0[,P1]
          # posterior hyperparameters
          dfN <- N-K-P+1
          ScaleN <- S[P1,P1]*tXXi/(N-K-P+1) # off-diagonal elements have no meaning
          meanN <- as.matrix(c(BetaHat[,P1]))
          # exclude inactive rows
          if(is.null(RrE[[h]])){RrE_h=NULL
          }else{
            if(nrow(RrE[[h]])==1){
              RrE_h <- t(as.matrix(RrE[[h]][,c((P1-1)*K+1:K,P*K+1)]))
            }else{
              RrE_h <- RrE[[h]][,c((P1-1)*K+1:K,P*K+1)]
            }
          }
          if(is.null(RrO[[h]])){RrO_h=NULL
          }else{
            if(nrow(RrO[[h]])==1){
              RrO_h <- t(as.matrix(RrO[[h]][,c((P1-1)*K+1:K,P*K+1)]))
            }else{
              RrO_h <- RrO[[h]][,c((P1-1)*K+1:K,P*K+1)]
            }
          }
          # compute relative measures of fit and complexity
          relcomp_h <- Student_measures(mean0,Scale0,df0,RrE_h,RrO_h,names1=matrixnames[,P1],
                                        constraints1=parse_hyp$original_hypothesis[h])
          relfit_h <- Student_measures(meanN,ScaleN,dfN,RrE_h,RrO_h,names1=matrixnames[,P1],
                                       constraints1=parse_hyp$original_hypothesis[h])

        }else{ #use Matrix-Student distributions with Monte Carlo estimate
          df0 <- 1
          dfN <- N-K-P+1
          relfit_h <- MatrixStudent_measures(Mean1=BetaHat,Scale1=S,tXXi1=tXXi,df1=dfN,RrE1=RrE[[h]],RrO1=RrO[[h]],
                                             Names1=matrix(names_coef,ncol=P),constraints1=parse_hyp$original_hypothesis[h],
                                             MCdraws=1e4)
          relcomp_h <- MatrixStudent_measures(Mean1=Mean0,Scale1=S_b,tXXi1=tXXi_b,df1=df0,RrE1=RrE[[h]],RrO1=RrO[[h]],
                                              Names1=matrix(names_coef,ncol=P),constraints1=parse_hyp$original_hypothesis[h],
                                              MCdraws=1e4)
        }
        return(list(relfit_h,relcomp_h))
      }))

      relfit <- t(matrix(unlist(relmeasunlist)[rep((0:(numhyp-1))*4,each=2)+rep(1:2,numhyp)],nrow=2))
      row.names(relfit) <- parse_hyp$original_hypothesis
      colnames(relfit) <- c("f_E","f_O")
      relcomp <- t(matrix(unlist(relmeasunlist)[rep((0:(numhyp-1))*4,each=2)+rep(3:4,numhyp)],nrow=2))
      row.names(relcomp) <- parse_hyp$original_hypothesis
      colnames(relcomp) <- c("c_E","c_O")

      # Compute relative fit/complexity for the complement hypothesis
      relfit <- MatrixStudent_prob_Hc(BetaHat,S,tXXi,N-K-P+1,as.matrix(relfit),RrO)
      relcomp <- MatrixStudent_prob_Hc(Mean0,S_b,tXXi_b,1,as.matrix(relcomp),RrO)
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
      row.names(BFtable) <- names(PHP_confirmatory)
      colnames(BFtable) <- c("comp_E","comp_O","fit_E","fit_O","BF_E","BF_O","BF","PHP")
      BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory))/
        t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
      row.names(BFmatrix_confirmatory) <- colnames(BFmatrix_confirmatory) <- names(BFtu_confirmatory)
      if(nrow(relfit)==length(parse_hyp$original_hypothesis)){
        hypotheses <- parse_hyp$original_hypothesis
      }else{
        hypotheses <- c(parse_hyp$original_hypothesis,"complement")
      }

    }

  }else{
    BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relfit <-
      relcomp <- hypotheses <- BFtable <- priorprobs <- NULL
  }

  bayesfactor <- "generalized adjusted fractional Bayes factor for coefficients;
    Bayes factor based on joint uniform priors for correlations"
  testedparameter <- "regression coefficients; correlations"

  BFlm_out <- list(
    BFtu_exploratory=BFtu_exploratory,
    PHP_exploratory=PHP_exploratory,
    BFtu_confirmatory=BFtu_confirmatory,
    PHP_confirmatory=PHP_confirmatory,
    BFmatrix_confirmatory=BFmatrix_confirmatory,
    BFtable_confirmatory=BFtable,
    BFtu_main=BFtu_main,
    PHP_main=PHP_main,
    BFtu_interaction=PHP_interaction,
    PHP_interaction=PHP_interaction,
    prior=priorprobs,
    hypotheses=hypotheses,
    estimates=postestimates,
    model=x,
    bayesfactor=bayesfactor,
    parameter=testedparameter,
    call=match.call())

  class(BFlm_out) <- "BF"

  return(BFlm_out)
}



params_in_hyp <- function(hyp){
  params_in_hyp <- trimws(unique(strsplit(hyp, split = "[ =<>,\\(\\);&\\*+-]+", perl = TRUE)[[1]]))
  params_in_hyp <- params_in_hyp[!sapply(params_in_hyp, grepl, pattern = "^[0-9]*\\.?[0-9]+$")]
  params_in_hyp[grepl("^[a-zA-Z]", params_in_hyp)]
}

#dyn.load("/Users/jorismulder/surfdrive/R packages/BFpack/scr/bct_continuous_final.dll")
# R function to call Fortran subroutine for Gibbs sampling using noninformative improper
# priors for regression coefficients, Jeffreys priors for standard deviations, and a proper
# joint uniform prior for the correlation matrices.
# dyn.load("/Users/jorismulder/surfdrive/R packages/BFpack/src/bct_continuous_final.dll")
estimate_postMeanCov_FisherZ <- function(YXlist,numdraws=5e3){
  # YXlist should be a list of length number of independent groups, of which each
  # element is another list of which the first element is a matrix of dependent
  # variables in the group, and the second element is a matrix of covariate variables
  # in the group.

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
  samsize0 <- numdraws
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
                 seed=as.integer(sample.int(1e6, 1)))

  FmeansCovCorr <- lapply(1:numG,function(g){
    Fdraws_g <- FisherZ(t(matrix(unlist(lapply(1:samsize0,function(s){
      res$CDrawsStore[s,g,,][lower.tri(diag(P))]
    })),ncol=samsize0)))
    mean_g <- apply(Fdraws_g,2,mean)
    covm_g <- cov(Fdraws_g)
    return(list(mean_g,covm_g))
  })
  meansCovCorr <- lapply(1:numG,function(g){
    draws_g <- t(matrix(unlist(lapply(1:samsize0,function(s){
      res$CDrawsStore[s,g,,][lower.tri(diag(P))]
    })),ncol=samsize0))
    mean_g <- apply(draws_g,2,mean)
    return(mean_g)
  })
  meanN <- unlist(lapply(1:numG,function(g){
    FmeansCovCorr[[g]][[1]]
  }))
  covmN <- matrix(0,nrow=numcorr,ncol=numcorr)
  numcorrg <- numcorr/numG
  for(g in 1:numG){
    covmN[(g-1)*numcorrg+1:numcorrg,(g-1)*numcorrg+1:numcorrg] <- FmeansCovCorr[[g]][[2]]
  }
  return(list(corr_quantiles=res$C_quantiles,B_quantiles=res$B_quantiles,
              sigma_quantiles=res$sigma_quantiles,meanN=meanN,covmN=covmN,
              corr_means=meansCovCorr))
}


# Fisher Z tranformation for sampled correlations
FisherZ <- function(r){.5*log((1+r)/(1-r))}


