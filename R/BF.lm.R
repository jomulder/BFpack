#' @importFrom pracma rref Rank
#' @importFrom mvtnorm dmvnorm pmvnorm rmvnorm dmvt pmvt rmvt
#' @importFrom stats rWishart qt
#' @importFrom MASS ginv
#' @describeIn BF S3 method for an object of class 'lm'
#' @method BF lm
#' @export
BF.lm <- function(x,
                  hypothesis = NULL,
                  prior.hyp = NULL,
                  complement = TRUE,
                  log = FALSE,
                  BF.type = 2,
                  ...){

  if(is.null(BF.type)){
    stop("The argument 'BF.type' must be the integer 1 (for the fractional BF) or 2 (for the adjusted fractional BF).")
  }
  if(!is.null(BF.type)){
    if(is.na(BF.type) | (BF.type!=1 & BF.type!=2))
      stop("The argument 'BF.type' must be the integer 1 (for the fractional BF) or 2 (for the adjusted fractional BF).")
  }
  if(BF.type==2){
    bayesfactor <- "generalized adjusted fractional Bayes factors"
  }else{
    bayesfactor <- "generalized fractional Bayes factors"
  }
  testedparameter <- "regression coefficients"

  logIN <- log

  # default BF on location parameters in a univarite normal linear model
  # Note that it is recommended that the fitten model is based on standardized covariates.
  x$coefficients <- as.matrix(x$coefficients)
  x$residuals <- as.matrix(x$residuals)

  P <- ncol(x$residuals)
  N <- nrow(x$residuals)
  K <- length(x$coefficients)/P # dimension of predictors per dependent variable
  dummyX <- rep(F,K)
  names(dummyX) <- row.names(x$coefficients)

  Xmat <- model.matrix(x)
  Ymat <- model.matrix(x)%*%x$coefficients + x$residuals

  # Exploratory testing of regression coefficients
  if(length(x$xlevels)==0){ #no grouping covariates: all observations receive equal fractions
    J <- 1
    dummyX <- rep(F,K)
    Nj <- nrow(Xmat)
    dvec <- rep(1,Nj)
    #set minimal fractions for the group
    bj <- ((P+K)/J)/Nj
    #no dummy covariates for factors
    dummy01TRUE <- FALSE
  }else{
    # factors are present, which are used to identify groups, which are used to determine the fractions
    # fractions are chosen such get minimal information is used per unique group (which is defined by the
    # unique combination of group indicators)

    numlevels <- unlist(lapply(x$xlevels,length))
    main_names <- names(x$xlevels)
    names_coef <- colnames(model.matrix(x))

    dummyX <- apply(matrix(unlist(lapply(1:length(main_names),function(faclev){
      unlist(lapply(1:length(names_coef),function(cf){
        grepl(main_names[faclev],names_coef[cf],fixed=TRUE)
      }))
    })),nrow=length(names_coef)),1,max)==1
    #number of groups on variations of dummy combinations
    groupcode <- as.matrix(unique(Xmat[,dummyX]))

    N <- nrow(Xmat)

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
    if(max(bj)>1){#then too few observations in certain groups, then use one single minimal fraction
      bj <- rep((P+K)/sum(Nj),length=J)
      if(bj[1]>1){
        stop("Not enough observations to compute a fractional Bayes factor.")
      }
    }

    # dummy01TRUE <- prod(unlist(lapply(1:length(x$contrasts),function(fac){
    #   x$contrasts[[fac]] == "contr.treatment"
    # }))) == 1
    # if(dummy01TRUE){
    #   numlevels <- unlist(lapply(x$xlevels,length))
    #   mains <- unlist(lapply(1:length(x$xlevels),function(fac){
    #     unlist(lapply(1:length(x$xlevels[[fac]]),function(lev){
    #       paste0(names(x$xlevels)[fac],x$xlevels[[fac]][lev])
    #     }))
    #   }))
    #   intercept <- attr(x$terms,"intercept")==1
    #   names_coef <- row.names(x$coefficients)
    #   # dummyX1 checks which columns of the design matrix X are dummy's for a
    #   # main effect or interaction effect
    #   dummyX1 <- apply(matrix(unlist(lapply(1:length(mains),function(faclev){
    #     unlist(lapply(1:length(names_coef),function(cf){
    #       grepl(mains[faclev],names_coef[cf],fixed=TRUE)
    #     }))
    #   })),nrow=length(names_coef)),1,max)==1
    # }else{
    #   dummyX1 <- rep(TRUE,K)
    # }
    # # dummyX2 checks which columns of the design matrix have two possible outcomes,
    # # which indicates a dummy variable
    # dummyX2 <- unlist(lapply(1:K,function(k){
    #   length(table(Xmat[,k])) == 2
    # }))
    # # dummyX indicate which columns contain dummy group covariates
    # dummyX <- dummyX2 * dummyX1 == 1
    # #number of groups on variations of dummy combinations
    # groupcode <- as.matrix(unique(Xmat[,dummyX]))
    # rownames(groupcode) <- unlist(lapply(1:nrow(groupcode),function(r){
    #   paste0("groupcode",r)
    # }))
    # J <- nrow(groupcode)
    # if(J==nrow(Xmat)){
    #   stop("Not enough observations for every group. Try fitting the model without factors.")
    # }
    #
    # # group membership of each observation
    # dvec <- unlist(lapply(1:N,function(i){
    #   which(rowSums(abs(t(matrix(rep(Xmat[i,dummyX],J),ncol=J)) - groupcode))==0)
    # }))
    # Nj <- c(table(dvec))
    # #set minimal fractions for each group
    # bj <- ((P+K)/J)/Nj
    # if(max(bj)>1){#then too few observations in certain groups, then use one single minimal fraction
    #   bj <- rep((P+K)/sum(Nj),length=J)
    #   if(bj[1]>1){
    #     stop("Not enough observations to compute a fractional Bayes factor.")
    #   }
    # }
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
      as.matrix(Xmat[dvec==j,])%*%t(Ymat[dvec==j,])
    } else {t(Xmat[dvec==j,])%*%Ymat[dvec==j,]}
  })
  tXYj_b <- lapply(1:J,function(j){
    tXYj[[j]]*bj[j]
  })
  tYYj <- lapply(1:J,function(j){
    if(Nj[j]==1){
      Ymat[dvec==j,]%*%t(Ymat[dvec==j,])
    }else t(Ymat[dvec==j,])%*%Ymat[dvec==j,]

  })
  tYYj_b <- lapply(1:J,function(j){
    tYYj[[j]]*bj[j]
  })
  tXX <- Reduce("+",tXXj)
  if(min(eigen(tXX)$values)<0){
    stop("Model matrix does not seem to be of full row rank.")
  }
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
  if(P==1){
    names_coef <- row.names(x$coefficients)
  }else{
    names_coef1 <- row.names(x$coefficients)
    names_coef2 <- colnames(x$coefficients)
    names_coef <- unlist(lapply(1:P,function(p){
      lapply(1:K,function(k){
        paste0(names_coef1[k],"_on_",names_coef2[p])
      })
    }))
  }

  # prior hyperparameters
  df0 <- 1 # should be the same as sum(rep(bj,times=Nj))-K-P+1
  Scale0 <- kronecker(S_b,tXXi_b)
  if(BF.type==2){
    mean0 <- as.matrix(rep(0,K*P))
  }else{
    mean0 <- as.matrix(c(BetaHat_b))
  }
  # posterior hyperparameters
  dfN <- N-K-P+1
  ScaleN <- kronecker(S,tXXi)/(N-K-P+1) # off-diagonal elements have no meaning
  meanN <- as.matrix(c(BetaHat))
  row.names(meanN) <- row.names(mean0) <- names_coef

  # Hypotheses for exploratory test
  # H0: beta = 0
  # H1: beta < 0
  # H2: beta > 0
  relfit <- t(matrix(unlist(lapply(1:(K*P),function(k){
    c(dt((0-meanN[k,1])/sqrt(ScaleN[k,k]),df=dfN,log=TRUE) - log(sqrt(ScaleN[k,k])),
      pt((0-meanN[k,1])/sqrt(ScaleN[k,k]),df=dfN,lower.tail = TRUE,log.p=TRUE),
      pt((0-meanN[k,1])/sqrt(ScaleN[k,k]),df=dfN,lower.tail = FALSE,log.p=TRUE))
  })),nrow=3))
  relcomp <- t(matrix(unlist(lapply(1:(K*P),function(k){
    c(dt((0-mean0[k,1])/sqrt(Scale0[k,k]),df=df0,log=TRUE) - log(sqrt(Scale0[k,k])),
      pt((0-mean0[k,1])/sqrt(Scale0[k,k]),df=df0,lower.tail = TRUE,log.p=TRUE),
      pt((0-mean0[k,1])/sqrt(Scale0[k,k]),df=df0,lower.tail = FALSE,log.p=TRUE))
  })),nrow=3))
  colnames(relfit) <- colnames(relcomp) <- c("p(=0)","Pr(<0)","Pr(>0)")
  row.names(relcomp) <- row.names(relfit) <- names_coef

  BFtu_exploratory <- relfit - relcomp
  colnames(BFtu_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
  maxrows <- apply(BFtu_exploratory,1,max)
  PHP_exploratory <- exp(BFtu_exploratory - maxrows %*% t(rep(1,3))) /
    apply(exp(BFtu_exploratory - maxrows %*% t(rep(1,3))),1,sum)

  #compute estimates
  postestimates <- cbind(meanN,meanN,
                         t(matrix(unlist(lapply(1:length(meanN),function(coef){
                           ub <- qt(p=.975,df=dfN)*sqrt(ScaleN[coef,coef])+meanN[coef,1]
                           lb <- qt(p=.025,df=dfN)*sqrt(ScaleN[coef,coef])+meanN[coef,1]
                           return(c(ub,lb))
                         })),nrow=2))
  )
  row.names(postestimates) <- names_coef
  colnames(postestimates) <- c("mean","median","2.5%","97.5%")

  # Additional exploratory tests of main effects and interaction effects
  # in the case of an aov type object
  if(sum(class(x)=="aov")==1 & J > 1){
    testedparameter <- "group means"

    mains <- unlist(lapply(1:length(x$xlevels),function(fac){
      unlist(lapply(1:length(x$xlevels[[fac]]),function(lev){
        paste0(names(x$xlevels)[fac],x$xlevels[[fac]][lev])
      }))
    }))

    # check main effects
    BFmain <- unlist(lapply(1:length(numlevels),function(fac){
      # if( number of levels of factor is less than the number of columns of main effect in model.matrix
      # then to test the main effect all effects of these columns should be zero, and if the number of
      # levels is equal to the number of columns then the effects of the columns should be equal.)
      name1 <- names(numlevels[fac])
      mains1 <- mains[sum(numlevels[1:fac])-numlevels[fac]+1:numlevels[fac]]
      which0 <- unlist(lapply(1:length(colnames(Xmat)),function(col){
        sum(colnames(Xmat)[col]==mains1)==1
      }))
      if(sum(which0) > 0){
        if(numlevels[fac]==sum(which0)){
          #there is no reference group
          #test whether the effects of all dummys are equal
          if(P > 1){
            RrE_f <- matrix(0,nrow=sum(which0)-1,ncol=length(colnames(Xmat)))
            for(r1 in 1:(sum(which0)-1)){
              RrE_f[r1,which(which0)[r1]] <- 1
              RrE_f[r1,which(which0)[r1+1]] <- -1
            }
            RrE_f <- cbind(kronecker(diag(P),RrE_f),0)
            relcomp_f <- MatrixStudent_measures(Mean1=matrix(mean0,ncol=P),Scale1=S_b,tXXi1=tXXi_b,
                                                df1=df0,RrE1=RrE_f,RrO1=NULL,Names1=NULL,constraints1=NULL,
                                                MCdraws=1e4)
            relfit_f <- MatrixStudent_measures(Mean1=BetaHat,Scale1=S,tXXi1=tXXi,df1=dfN,RrE1=RrE_f,
                                               RrO1=NULL,Names1=NULL,constraints1=NULL,MCdraws=1e4)
          }else{
            RrE_f <- matrix(0,nrow=sum(which0)-1,ncol=length(colnames(Xmat))+1)
            for(r1 in 1:(sum(which0)-1)){
              RrE_f[r1,which(which0)[r1]] <- 1
              RrE_f[r1,which(which0)[r1+1]] <- -1
            }
            relcomp_f <- Student_measures(mean1=mean0,Scale1=Scale0,df1=df0,RrE1=RrE_f,RrO1=NULL)
            relfit_f <- Student_measures(mean1=meanN,Scale1=ScaleN,df1=dfN,RrE1=RrE_f,RrO1=NULL)
          }

        }else{
          #there is a reference group
          #test whether the effects of the remaining dummys are zero of these indicators equals zero
          if(P > 1){
            RrE_f <- matrix(0,nrow=sum(which0),ncol=length(colnames(Xmat)))
            for(r1 in 1:sum(which0)){
              RrE_f[r1,which(which0)[r1]] <- 1
            }
            RrE_f <- cbind(kronecker(diag(P),RrE_f),rep(0,sum(which0)*P))
            relcomp_f <- MatrixStudent_measures(Mean1=matrix(mean0,ncol=P),Scale1=S_b,tXXi1=tXXi_b,
                                                df1=df0,RrE1=RrE_f,RrO1=NULL,Names1=NULL,constraints1=NULL,
                                                MCdraws=1e4)
            relfit_f <- MatrixStudent_measures(Mean1=BetaHat,Scale1=S,tXXi1=tXXi,df1=dfN,RrE1=RrE_f,
                                               RrO1=NULL,Names1=NULL,constraints1=NULL,MCdraws=1e4)
          }else{
            RrE_f <- matrix(0,nrow=sum(which0),ncol=length(colnames(Xmat))+1)
            for(r1 in 1:sum(which0)){
              RrE_f[r1,which(which0)[r1]] <- 1
            }
            relcomp_f <- Student_measures(mean1=mean0,Scale1=Scale0,df1=df0,RrE1=RrE_f,RrO1=NULL)
            relfit_f <- Student_measures(mean1=meanN,Scale1=ScaleN,df1=dfN,RrE1=RrE_f,RrO1=NULL)
          }
        }
        BFtu <- relfit_f[1] - relcomp_f[1]
        names(BFtu) <- name1
        return(c(BFtu,relfit_f[1],relcomp_f[1]))
      }
    }))
    #compute Bayes factors for testing main effects if present
    if(length(BFmain)>0){ # then there are main effects
      names_main <- names(BFmain[(0:(length(BFmain)/3-1))*3+1])
      BFtu_main <- matrix(c(BFmain[(0:(length(BFmain)/3-1))*3+1],rep(log(1),length(BFmain)/3)),
                          nrow=length(BFmain)/3)
      row.names(BFtu_main) <- names_main
      colnames(BFtu_main) <- c("BFtu","BFuu")
      maxBFtu <- apply(BFtu_main,1,max)
      PHP_main <- exp(BFtu_main - maxBFtu %*% t(rep(1,ncol(BFtu_main)))) /
        apply(exp(BFtu_main - maxBFtu %*% t(rep(1,ncol(BFtu_main)))),1,sum)
      colnames(PHP_main) <- c("Pr(no effect)","Pr(full model)")
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

    for(c1 in 1:ncol(matcov)){
      if(c1 < ncol(matcov)){
        numeffects_c <- sum(matcov[,c1])
        if(numeffects_c>1){
          count_interaction <- count_interaction + 1
          interactionset <- names(which(matcov[,c1]))
          whichx <- apply(matcov[which(matcov[,c1]),],2,sum)==length(interactionset)
          if(P > 1){
            RrE_ia <- matrix(0,nrow=sum(whichx),ncol=K)
            for(r1 in 1:sum(whichx)){RrE_ia[r1,which(whichx)[r1]]<-1}
            RrE_ia <- cbind(kronecker(diag(P),RrE_ia),rep(0,sum(whichx)*P))
            relcomp_ia <- MatrixStudent_measures(Mean1=matrix(mean0,ncol=P),Scale1=S_b,tXXi1=tXXi_b,
                                                 df1=df0,RrE1=RrE_ia,RrO1=NULL,Names1=NULL,
                                                 constraints1=NULL,MCdraws=1e4)
            relfit_ia <- MatrixStudent_measures(Mean1=BetaHat,Scale1=S,tXXi1=tXXi,df1=dfN,RrE1=RrE_ia,
                                                RrO1=NULL,Names1=NULL,constraints1=NULL,MCdraws=1e4)
          }else{
            RrE_ia <- matrix(0,nrow=sum(whichx),ncol=K+1)
            for(r1 in 1:sum(whichx)){RrE_ia[r1,which(whichx)[r1]]<-1}
            relcomp_ia <- Student_measures(mean1=mean0,Scale1=Scale0,df1=df0,RrE1=RrE_ia,RrO1=NULL)
            relfit_ia <- Student_measures(mean1=meanN,Scale1=ScaleN,df1=dfN,RrE1=RrE_ia,RrO1=NULL)
          }
          names(relcomp_ia) <- c("c=","c>")
          names(relfit_ia) <- c("f=","f>")
          BFtu_ia <- relfit_ia[1] - relcomp_ia[1]
          names(BFtu_ia) <- paste(interactionset,collapse=":")
          BFtu_interaction0[[count_interaction]] <- c(BFtu_ia,relcomp_ia[1],relfit_ia[1])
          #exclude other columns from Xmat that have been used to avoid double results
          matcov <- matcov[,-(which(whichx[-c1])+1)]
        }
      }
    }
    #compute Bayes factors for testing interaction effects if present
    if(count_interaction>0){ # then there are interaction effects
      BFtu_interaction0 <- unlist(BFtu_interaction0)
      names_interaction <- names(BFtu_interaction0[(0:(length(BFtu_interaction0)/3-1))*3+1])
      BFtu_interaction <- matrix(c(BFtu_interaction0[(0:(length(BFtu_interaction0)/3-1))*3+1],
                                   rep(log(1),length(BFtu_interaction0)/3)),nrow=length(BFtu_interaction0)/3)
      row.names(BFtu_interaction) <- names_interaction
      colnames(BFtu_interaction) <- c("BFtu","BFuu")
      maxrows <- apply(BFtu_interaction,1,max)
      PHP_interaction <- exp(BFtu_interaction - maxrows %*% t(rep(1,ncol(BFtu_interaction)))) /
        apply(exp(BFtu_interaction - maxrows %*% t(rep(1,ncol(BFtu_interaction)))),1,sum)
      colnames(PHP_interaction) <- c("Pr(no effect)","Pr(full model)")
    }else{ PHP_interaction <- BFtu_interaction <- NULL}
    #BFtu_exploratory <- rbind(BFtu_main,BFtu_interaction)
    #PHP_exploratory <- rbind(PHP_main,PHP_interaction)
  }

  if(logIN == FALSE){
    BFtu_exploratory <- exp(BFtu_exploratory)
  }

  # confirmatory BF test
  if(!is.null(hypothesis)){

    #then constraints on regression coefficients
    matrixnames <- matrix(names_coef,nrow=K)

    # translate named constraints to matrices with coefficients for constraints
    parse_hyp <- bain:::parse_hypothesis(names_coef,hypothesis)
    parse_hyp$hyp_mat <- do.call(rbind, parse_hyp$hyp_mat)
    RrList <- make_RrList2(parse_hyp)
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]

    if(length(RrE)==1){
      RrStack <- rbind(do.call(rbind,RrE),do.call(rbind,RrO))
      RrStack <- interval_RrStack(RrStack)
    }else{
      RrStack_list <- lapply(1:length(RrE),function(h){
        interval_RrStack(rbind(RrE[[h]],RrO[[h]]))
      })
      RrStack <- do.call(rbind,RrStack_list)
    }
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
      nonzero <- rref_ei[,P*K+1]!=0
      if(max(nonzero)>0){
        row1 <- max(which(nonzero))
        if(sum(abs(rref_ei[row1,1:(P*K)]))==0){
          stop("No common boundary point for prior location. Conflicting constraints.")
        }
      }
    }

    #number of hypotheses that are specified
    numhyp <- length(RrO)

    if(P > 1){

      #default prior location
      if(BF.type==2){
        Mean0 <- matrix(c(ginv(RStack)%*%rStack),nrow=K,ncol=P)
      }else{
        Mean0 <- BetaHat #mean0
      }
      #Mean0 <- matrix(c(ginv(RStack)%*%rStack),nrow=K,ncol=P)

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
      if(complement==TRUE){
        relfit <- MatrixStudent_prob_Hc(Mean1=BetaHat,Scale1=S,tXXi1=tXXi,df1=N-K-P+1,
                                        relmeas=as.matrix(relfit),RrO1=RrO)
        relcomp <- MatrixStudent_prob_Hc(Mean0,S_b,tXXi_b,1,as.matrix(relcomp),RrO)
        hypothesisshort <- unlist(lapply(1:nrow(relfit),function(h) paste0("H",as.character(h))))
        row.names(relcomp) <- row.names(relfit) <- hypothesisshort
      }

      # # the BF for the complement hypothesis vs Hu needs to be computed.
      # BFtu_conf <- relfit - relcomp
      # BFtu_confirmatory <- c(apply(BFtu_conf, 1, sum))
      # # Check input of prior probabilies
      # if(is.null(prior.hyp)){
      #   priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
      # }else{
      #   if(!is.numeric(prior.hyp) || length(prior.hyp)!=length(BFtu_confirmatory)){
      #     warning(paste0("Argument 'prior.hyp' should be numeric and of length ",as.character(length(BFtu_confirmatory)),". Equal prior probabilities are used."))
      #     priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
      #   }else{
      #     priorprobs <- prior.hyp
      #   }
      # }
      #
      # names(priorprobs) <- names(BFtu_confirmatory)
      # BFtu_scaled <- BFtu_confirmatory - max(BFtu_confirmatory)
      # PHP_confirmatory <- exp(BFtu_scaled)*priorprobs / sum(exp(BFtu_scaled)*priorprobs)
      # BFtable <- cbind(relcomp,relfit,relfit[,1]/relcomp[,1],relfit[,2]/relcomp[,2],
      #                  apply(relfit,1,prod)/apply(relcomp,1,prod),PHP_confirmatory)
      # row.names(BFtable) <- names(PHP_confirmatory)
      # colnames(BFtable) <- c("comp_E","comp_O","fit_E","fit_O","BF_E","BF_O","BF","PHP")
      # BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)) -
      #   t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
      # diag(BFmatrix_confirmatory) <- 1
      # row.names(BFmatrix_confirmatory) <- colnames(BFmatrix_confirmatory) <- names(BFtu_confirmatory)
      # if(nrow(relfit)==length(parse_hyp$original_hypothesis)){
      #   hypotheses <- parse_hyp$original_hypothesis
      # }else{
      #   hypotheses <- c(parse_hyp$original_hypothesis,"complement")
      # }

    }else{
      # one dependent variable and posterior/prior have Student t distributions

      # prior hyperparameters
      df0 <- 1 # should be the same as sum(rep(bj,times=Nj))-K-P+1
      Scale0 <- kronecker(S_b,tXXi_b)
      #default prior location
      if(BF.type==2){
        mean0 <- ginv(RStack)%*%rStack
      }else{
        mean0 <- BetaHat #mean0
      }
      # posterior hyperparameters
      dfN <- N-K-P+1
      ScaleN <- kronecker(S,tXXi)/(N-K-P+1) # off-diagonal elements have no meaning
      meanN <- as.matrix(c(BetaHat))

      relcomp <- t(matrix(unlist(lapply(1:numhyp,function(h){
        Student_measures(mean1=mean0,Scale1=Scale0,df1=df0,RrE1=RrE[[h]],RrO1=RrO[[h]],
                         names1=names_coef,constraints1=parse_hyp$original_hypothesis[h])
      })),nrow=2))

      relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
        Student_measures(meanN,ScaleN,dfN,RrE[[h]],RrO[[h]],
                         names1=names_coef,constraints1=parse_hyp$original_hypothesis[h])
      })),nrow=2))
      colnames(relcomp) <- c("c_E","c_O")
      colnames(relfit) <- c("f_E","f_O")
      row.names(relcomp)[1:numhyp] <- parse_hyp$original_hypothesis
      row.names(relfit)[1:numhyp] <- parse_hyp$original_hypothesis

      if(complement == TRUE){
        # Compute relative fit/complexity for the complement hypothesis
        relfit <- Student_prob_Hc(meanN,ScaleN,dfN,relfit,hypothesis,RrO)
        relcomp <- Student_prob_Hc(mean0,Scale0,df0,relcomp,hypothesis,RrO)
        row.names(relcomp)[1:numhyp] <- parse_hyp$original_hypothesis
        row.names(relfit)[1:numhyp] <- parse_hyp$original_hypothesis
      }

      # # the BF for the complement hypothesis vs Hu needs to be computed.
      # BFtu_confirmatory <- c(apply(relfit - relcomp, 1, sum))
      # # Check input of prior probabilies
      # if(is.null(prior.hyp)){
      #   priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
      # }else{
      #   if(!is.numeric(prior.hyp) || length(prior.hyp)!=length(BFtu_confirmatory)){
      #     warning(paste0("Argument 'prior.hyp' should be numeric and of length ",as.character(length(BFtu_confirmatory)),". Equal prior probabilities are used."))
      #     priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
      #   }else{
      #     priorprobs <- prior.hyp
      #   }
      # }
      # names(priorprobs) <- names(BFtu_confirmatory)
      # PHP_confirmatory <- BFtu_confirmatory*priorprobs / sum(BFtu_confirmatory*priorprobs)
      # BFtable <- cbind(relcomp,relfit,relfit[,1]/relcomp[,1],relfit[,2]/relcomp[,2],
      #                  BFtu_confirmatory,PHP_confirmatory)
      # row.names(BFtable) <- names(BFtu_confirmatory)
      # colnames(BFtable) <- c("complex=","complex>","fit=","fit>","BF=","BF>","BF","PHP")
      # BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory))/
      #   t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
      # diag(BFmatrix_confirmatory) <- 1
      # row.names(BFmatrix_confirmatory) <- colnames(BFmatrix_confirmatory) <- names(BFtu_confirmatory)
      # #tested hypotheses
      # hypotheses <- row.names(relfit)
    }

    # the BF for the complement hypothesis vs Hu needs to be computed.
    BFtu_conf <- relfit - relcomp
    BFtu_confirmatory <- c(apply(BFtu_conf, 1, sum))
    # Check input of prior probabilies
    if(is.null(prior.hyp)){
      priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
    }else{
      if(!is.numeric(prior.hyp) || length(prior.hyp)!=length(BFtu_confirmatory)){
        warning(paste0("Argument 'prior.hyp' should be numeric and of length ",as.character(length(BFtu_confirmatory)),". Equal prior probabilities are used."))
        priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
      }else{
        priorprobs <- prior.hyp
      }
    }

    names(priorprobs) <- names(BFtu_confirmatory)
    BFtu_scaled <- BFtu_confirmatory - max(BFtu_confirmatory)
    PHP_confirmatory <- exp(BFtu_scaled)*priorprobs / sum(exp(BFtu_scaled)*priorprobs)
    BFtable <- cbind(relcomp,relfit,relfit[,1]-relcomp[,1],relfit[,2]-relcomp[,2],
                     apply(relfit,1,sum)-apply(relcomp,1,sum),PHP_confirmatory)
    BFtable[,1:7] <- exp(BFtable[,1:7])
    row.names(BFtable) <- names(PHP_confirmatory)
    colnames(BFtable) <- c("comp_E","comp_O","fit_E","fit_O","BF_E","BF_O","BF","PHP")
    BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)) -
      t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
    diag(BFmatrix_confirmatory) <- log(1)
    row.names(BFmatrix_confirmatory) <- colnames(BFmatrix_confirmatory) <- names(BFtu_confirmatory)
    if(nrow(relfit)==length(parse_hyp$original_hypothesis)){
      hypotheses <- parse_hyp$original_hypothesis
    }else{
      hypotheses <- c(parse_hyp$original_hypothesis,"complement")
    }
    if(logIN == FALSE){
      BFtu_confirmatory <- exp(BFtu_confirmatory)
      BFmatrix_confirmatory <- exp(BFmatrix_confirmatory)
    }

  }else{
    BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relfit <-
      relcomp <- hypotheses <- BFtable <- priorprobs <- NULL
  }

  if(sum(class(x)=="aov")==1 & J > 1){
    BFlm_out <- list(
      BFtu_exploratory=BFtu_exploratory,
      BFtu_main=BFtu_main,
      BFtu_interaction=BFtu_interaction,
      PHP_exploratory=PHP_exploratory,
      PHP_main=PHP_main,
      PHP_interaction=PHP_interaction,
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
      log = logIN,
      fraction_number_groupIDs = J,
      fraction_groupID_observations = dvec,
      call=match.call())
  }else{
    BFlm_out <- list(
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
      log = logIN,
      fraction_number_groupIDs = J,
      fraction_groupID_observations = dvec,
      call=match.call())
  }

  names(BFlm_out$fraction_group_identifier) <- NULL

  class(BFlm_out) <- "BF"

  return(BFlm_out)
}



# compute relative meausures (fit or complexity) under a multivariate Student t distribution
MatrixStudent_measures <- function(Mean1,Scale1,tXXi1,df1,RrE1,RrO1,Names1=NULL,
                                   constraints1=NULL,MCdraws=1e4){
  # constraints1 = parse_hyp$original_hypothesis
  # RrE1 <- matrix(0,nrow=1,ncol=ncol(RrO1))
  # RrE1[1,1:2] <- c(-1,1)
  K <- nrow(Mean1)
  P <- ncol(Mean1)
  # vectorize the mean
  mean1 <- c(Mean1)
  relE <- relO <- log(1)
  if(!is.null(RrE1) && is.null(RrO1)){ #only equality constraints
    RE1 <- RrE1[,-(K*P+1)]
    if(!is.matrix(RE1)){
      RE1 <- matrix(RE1,ncol=K*P)
    }
    rE1 <- RrE1[,(K*P+1)]
    qE1 <- nrow(RE1)

    temp1 <- rWishart(MCdraws,df1+P-1,solve(Scale1))
    temp2 <- lapply(seq(dim(temp1)[3]), function(x) temp1[,,x])
    SigmaList <- lapply(temp2,solve)
    covm1_E <- lapply(SigmaList,function(temp) RE1%*%(kronecker(temp,tXXi1))%*%t(RE1) )
    mean1_E <- c(RE1 %*% mean1)
    relE <- log(mean(unlist(lapply(covm1_E,function(temp) dmvnorm(rE1,mean=mean1_E,sigma=temp)))))

  }else{
    if(is.null(RrE1) && !is.null(RrO1)){ #only order constraints
      RO1 <- RrO1[,-(K*P+1)]
      if(!is.matrix(RO1)){
        RO1 <- matrix(RO1,ncol=K*P)
      }
      qO1 <- nrow(RO1)
      rO1 <- RrO1[,(K*P+1)]

      if(Rank(RO1)==nrow(RO1)){ #RO1 is of full row rank. So use transformation.

        Scale1inv <- solve(Scale1)
        relO <- unlist(lapply(1:1e3,function(s){
          Sigma1 <- solve(rWishart(1,df=df1+P-1,Sigma=Scale1inv)[,,1])
          meanO <- c(RO1%*%mean1)
          covmO <- RO1%*%kronecker(Sigma1,tXXi1)%*%t(RO1)
          pmvnorm(lower=rO1,upper=Inf,mean=meanO,sigma=covmO)[1]
        }))
        relO <- log(mean(relO[relO!="NaN"]))
        if(relO>0){relO <- 0}

      }else{ #no linear transformation can be used; pmvt cannot be used. Use bain with a multivariate normal approximation
        #compute covariance matrix for multivariate normal distribution
        mean1 <- c(Mean1)
        names(mean1) <- c(Names1)
        if(df1>2){ #posterior measures
          covm1 <- kronecker(Scale1,tXXi1)/(df1-2)
          bain_res <- bain(x=mean1,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
          relO <- log(bain_res$fit[1,3])
        }else if(df1==2){ #posterior measures
          covm1 <- kronecker(Scale1,tXXi1)/(df1-1)
          bain_res <- bain(x=mean1,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
          relO <- log(bain_res$fit[1,3])
        }else{
          covm1 <- kronecker(Scale1,tXXi1) #for prior with df1==1, probability independent of common factor of scale1
          bain_res <- bain(x=mean1,hypothesis=constraints1,Sigma=covm1,n=df1) #n not used in computation
          relO <- log(bain_res$fit[1,4])
        }

        # bain1 <- bain::bain(mean1,Sigma1=covm1,RrE1,RrO1,n=10) # choice of n does not matter
        # extract posterior probability (Fit_eq) from bain-object)
        # warning("Check if this works now")
      }
    }else{ #hypothesis with equality and order constraints

      RE1 <- RrE1[,-(K*P+1)]
      if(!is.matrix(RE1)){
        RE1 <- matrix(RE1,ncol=K*P)
      }
      rE1 <- RrE1[,(K*P+1)]
      qE1 <- nrow(RE1)
      RO1 <- RrO1[,-(K*P+1)]
      if(!is.matrix(RO1)){
        RO1 <- matrix(RO1,ncol=K*P)
      }
      qO1 <- nrow(RO1)
      rO1 <- RrO1[,(K*P+1)]
      Rr1 <- rbind(RrE1,RrO1)
      R1 <- rbind(RE1,RO1)
      r1 <- c(rE1,rO1)
      qC1 <- length(r1)

      temp1 <- rWishart(MCdraws,df1+P-1,solve(Scale1))
      temp2 <- lapply(seq(dim(temp1)[3]), function(x) temp1[,,x])
      SigmaList <- lapply(temp2,solve)
      covm1_E <- lapply(SigmaList,function(temp) RE1%*%(kronecker(temp,tXXi1))%*%t(RE1) )
      mean1_E <- RE1 %*% mean1
      relE <- log(mean(unlist(lapply(covm1_E,function(temp) dmvnorm(rE1,mean=mean1_E,sigma=temp)))))

      if(Rank(Rr1) == nrow(Rr1)){
        covm1_O <- lapply(SigmaList,function(temp) R1%*%(kronecker(temp,tXXi1))%*%t(R1) )
        mean1_O <- c(R1%*%mean1)

        mean1_OE <- lapply(covm1_O,function(temp){
          as.vector(mean1_O[(qE1+1):qC1] +
                      c(matrix(temp[(qE1+1):qC1,1:qE1],ncol=qE1)%*%solve(temp[1:qE1,1:qE1])%*%(rE1-mean1_E)))
        })
        covm1_OE <- lapply(covm1_O,function(temp) temp[(qE1+1):qC1,(qE1+1):qC1] -
                             matrix(temp[(qE1+1):qC1,1:qE1],ncol=qE1)%*%solve(temp[1:qE1,1:qE1])%*%
                             matrix(temp[1:qE1,(qE1+1):qC1],nrow=qE1))
        #check covariance because some can be nonsymmetric due to a generation error
        welk1 <- which(unlist(lapply(covm1_OE,function(temp) isSymmetric(temp,
                                                                         tol = sqrt(.Machine$double.eps),check.attributes = FALSE) &&
                                       min(eigen(temp)$values)>sqrt(.Machine$double.eps) )))
        covm1_OE <- covm1_OE[welk1]
        mean1_OE <- mean1_OE[welk1]
        relO <- log(mean(mapply(function(mu_temp,Sigma_temp) pmvnorm(lower=rO1,
                                                                     upper=rep(Inf,qO1),mean=mu_temp,sigma=Sigma_temp)[1],mean1_OE,covm1_OE)))
        if(relO > 0){relO <- 0}
      }else{ #use bain for the computation of the probability

        mean1 <- c(Mean1)
        names(mean1) <- c(Names1)
        if(df1>2){ #posterior measures
          covm1 <- kronecker(Scale1,tXXi1)/(df1-2)
          bain_res <- bain(x=mean1,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
          relO <- log(bain_res$fit[1,3])
        }else if(df1==2){ #posterior measures
          covm1 <- kronecker(Scale1,tXXi1)/(df1-1)
          bain_res <- bain(x=mean1,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
          relO <- log(bain_res$fit[1,3])
        }else{
          covm1 <- kronecker(Scale1,tXXi1) #for prior with df1==1, probability independent of common factor of scale1
          bain_res <- bain(x=mean1,hypothesis=constraints1,Sigma=covm1,n=df1) #n not used in computation
          relO <- log(bain_res$fit[1,4])
        }
      }
    }
  }
  return(c(relE,relO))
}

# compute relative meausures (fit or complexity) under a multivariate Student t distribution
Student_measures <- function(mean1,Scale1,df1,RrE1,RrO1,names1=NULL,constraints1=NULL){ # Volgens mij moet je hier ook N meegeven
  K <- length(mean1)
  relE <- relO <- log(1)
  if(!is.null(RrE1) && is.null(RrO1)){ #only equality constraints
    RE1 <- RrE1[,-(K+1)]
    if(!is.matrix(RE1)){
      RE1 <- matrix(RE1,ncol=K)
    }
    rE1 <- RrE1[,(K+1)]
    qE1 <- nrow(RE1)
    meanE <- RE1%*%mean1
    scaleE <- RE1%*%Scale1%*%t(RE1)
    relE <- dmvt(rE1,delta=c(meanE),sigma=scaleE,df=df1,log=TRUE)
  }
  if(is.null(RrE1) && !is.null(RrO1)){ #only order constraints
    RO1 <- RrO1[,-(K+1)]
    if(!is.matrix(RO1)){
      RO1 <- matrix(RO1,ncol=K)
    }
    qO1 <- nrow(RO1)
    rO1 <- RrO1[,(K+1)]

    if(Rank(RO1)==nrow(RO1)){ #RO1 is of full row rank. So use transformation.
      meanO <- c(RO1%*%mean1)
      scaleO <- RO1%*%Scale1%*%t(RO1)
      if(nrow(scaleO)==1){
        relO <- pt((rO1-meanO)/sqrt(scaleO[1,1]),df=df1,lower.tail=FALSE,log.p = TRUE)
      }else{
        relO <- pmvt(lower=rO1,upper=Inf,delta=meanO,sigma=scaleO,df=df1,
             type="shifted",algorithm="TVPACK")[1]
        if(relO<0){relO <- 0}
        if(relO>1){relO <- 1}
        relO <- log(relO)
      }
    }else{ #no linear transformation can be used; pmvt cannot be used. Use bain with a multivariate normal approximation
      #compute covariance matrix for multivariate normal distribution
      row.names(mean1) <- names1
      if(df1>2){ # we need posterior measures
        covm1 <- Scale1/(df1-2)
        mean1vec <- c(mean1)
        names(mean1vec) <- row.names(mean1)
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
        relO <- log(bain_res$fit[1,3])
      }else if(df1==2){ # we need posterior measures (there is very little information)
        covm1 <- Scale1/(df1-1)
        mean1vec <- c(mean1)
        names(mean1vec) <- row.names(mean1)
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
        relO <- log(bain_res$fit[1,3])
      }else{ #then df=1, so we need prior measures
        covm1 <- Scale1 #for prior with df1==1, probability independent of common factor of scale1
        mean1vec <- c(mean1)
        names(mean1vec) <- row.names(mean1)
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=covm1,n=df1) #n not used in computation
        relO <- log(bain_res$fit[1,4])
      }
    }
  }
  if(!is.null(RrE1) && !is.null(RrO1)){ #hypothesis with equality and order constraints

    RE1 <- RrE1[,-(K+1)]
    if(!is.matrix(RE1)){
      RE1 <- matrix(RE1,ncol=K)
    }
    rE1 <- RrE1[,(K+1)]
    qE1 <- nrow(RE1)
    RO1 <- RrO1[,-(K+1)]
    if(!is.matrix(RO1)){
      RO1 <- matrix(RO1,ncol=K)
    }
    qO1 <- nrow(RO1)
    rO1 <- RrO1[,(K+1)]

    #a)Transformation matrix
    D <- diag(K) - t(RE1) %*% solve(RE1 %*% t(RE1)) %*% RE1
    D2 <- unique(round(D, 5))
    if(length(as.logical(rowSums(D2 != 0)))==1){
      D2 <- matrix(D2[as.logical(rowSums(D2 != 0)),],nrow=1)
    }else{
      D2 <- D2[as.logical(rowSums(D2 != 0)),]
    }
    if(!is.matrix(D2)){
      D2 <- t(D2)
    }
    Tm <- rbind(RE1, D2)

    #b)
    Tmean1 <- Tm %*% mean1
    Tscale1 <- Tm %*% Scale1 %*% t(Tm)

    # relative meausure for equalities
    relE <- dmvt(x = t(rE1), delta = Tmean1[1:qE1], sigma = matrix(Tscale1[1:qE1, 1:qE1],
                                                                   ncol = qE1), df = df1, log = TRUE)

    # transform order constraints
    RO1tilde <- RO1 %*% ginv(D2)
    rO1tilde <- rO1 - RO1 %*% ginv(RE1) %*% rE1

    # Partitioning equality part and order part
    Tmean1E <- Tmean1[1:qE1]
    Tmean1O <- Tmean1[(qE1 + 1):K]

    Tscale1EE <- as.matrix(Tscale1[1:qE1, 1:qE1],nrow=qE1)
    Tscale1OE <- matrix(Tscale1[(qE1 + 1):K, 1:qE1],ncol=qE1)
    Tscale1OO <- as.matrix(Tscale1[(qE1 + 1):K, (qE1 + 1):K],ncol=K-qE1)

    #conditional location and scale matrix
    Tmean1OgE <- Tmean1O + Tscale1OE %*% solve(Tscale1EE) %*% matrix(rE1 - Tmean1E)
    Tscale1OgE <- as.vector((df1 + (t(matrix(rE1 - Tmean1E)) %*% solve(Tscale1EE) %*% matrix(rE1 - Tmean1E))) /
                              (df1 + qE1)) * (Tscale1OO - Tscale1OE %*% solve(Tscale1EE) %*% t(Tscale1OE))

    if(Rank(RO1tilde) == nrow(RO1tilde)){
      rO1tilde <- as.vector(rO1tilde)

      delta_trans <- as.vector(RO1tilde %*% Tmean1OgE)
      scale1_trans <- RO1tilde %*% Tscale1OgE %*% t(RO1tilde)

      if(nrow(scale1_trans)==1){ # univariate
        relO <- pt((rO1tilde - delta_trans) / sqrt(scale1_trans), df = df1+qE1, lower.tail = FALSE,
                   log.p = TRUE)[1]
      } else { # multivariate
        relO <- log(pmvt(lower = rO1tilde, upper = Inf, delta = delta_trans, sigma = scale1_trans,
                         df = df1+qE1, type = "shifted",algorithm="TVPACK")[1])
        if(relO<0){relO <- 0}
        if(relO>1){relO <- 1}
        relO <- log(relO)
      }

    }else{ #use bain for the computation of the probability
      #compute covariance matrix for multivariate normal distribution
      row.names(mean1) <- names1
      if(df1>2){ # we need posterior measures
        covm1 <- Scale1/(df1-2)
        mean1vec <- c(mean1)
        names(mean1vec) <- row.names(mean1)
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
        relO <- log(bain_res$fit[1,3])
      }else if(df1==2){ # we need posterior measures (there is very little information)
        covm1 <- Scale1/(df1-1)
        mean1vec <- c(mean1)
        names(mean1vec) <- row.names(mean1)
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
        relO <- log(bain_res$fit[1,3])
      }else{ #then df=1, so we need prior measures
        covm1 <- Scale1 #for prior with df1==1, probability independent of common factor of scale1
        mean1vec <- c(mean1)
        names(mean1vec) <- row.names(mean1)
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=covm1,n=df1) #n not used in computation
        relO <- log(bain_res$fit[1,4])
      }
    }
  }

  return(c(relE,relO))
}

# The function computes the probability of an unconstrained draw falling in the complement subspace.
MatrixStudent_prob_Hc <- function(Mean1,Scale1,tXXi1,df1,relmeas,RrO1){

  P <- ncol(Mean1)
  K <- nrow(Mean1)
  numpara <- P*K
  numhyp <- nrow(relmeas)
  #  relmeas <- relmeas[1:numhyp,]
  which_eq <- relmeas[,1] != log(1)
  if(sum(which_eq)==numhyp){ # Then the complement is equivalent to the unconstrained hypothesis.
    relmeas <- rbind(relmeas,log(rep(1,2)))
    rownames(relmeas)[numhyp+1] <- "complement"
  }else{ # So there is at least one hypothesis with only order constraints
    welk <- which(!which_eq)
    if(length(welk)==1){ # There is one hypothesis with only order constraints. Hc is complement of this hypothesis.
      relmeas <- rbind(relmeas,rep(1,2))
      relmeas[numhyp+1,2] <- 1 - exp(relmeas[welk,2])
      relmeas[numhyp+1,] <- log(relmeas[numhyp+1,])
      rownames(relmeas)[numhyp+1] <- "complement"
    }else{ # So more than one hypothesis with only order constraints
      # First we check whether ther is an overlap between the order constrained spaces.

      draws2 <- 1e4
      randomDraws <- rmvnorm(draws2,mean=rep(0,numpara),sigma=diag(numpara))
      #get draws that satisfy the constraints of the separate order constrained hypotheses
      checksOC <- lapply(welk,function(h){
        Rorder <- as.matrix(RrO1[[h]][,-(1+numpara)])
        if(ncol(Rorder)==1){
          Rorder <- t(Rorder)
        }
        rorder <- as.matrix(RrO1[[h]][,1+numpara])
        apply(randomDraws%*%t(Rorder) > rep(1,draws2)%*%t(rorder),1,prod)
      })
      checkOCplus <- Reduce("+",checksOC)

      if(sum(checkOCplus > 0) < draws2){ #then the joint order constrained hypotheses do not completely cover the parameter space.
        if(sum(checkOCplus>1)==0){ # then order constrained spaces are nonoverlapping
          relmeas <- rbind(relmeas,rep(1,2))
          relmeas[numhyp+1,2] <- 1 - sum(exp(relmeas[welk,2]))
          relmeas[numhyp+1,] <- log(relmeas[numhyp+1,])
          rownames(relmeas)[numhyp+1] <- "complement"
        }else{ #the order constrained subspaces at least partly overlap

          # function below gives a rough estimate of the posterior probability under Hc
          # a bain type of algorithm would be better of course. but for now this is ok.

          temp1 <- rWishart(draws2,df1+P-1,solve(Scale1))
          temp2 <- lapply(seq(dim(temp1)[3]), function(x) temp1[,,x])
          SigmaList <- lapply(temp2,solve)
          randomDraws <- matrix(unlist(lapply(SigmaList,function(temp){
            rmvnorm(1,mean=c(Mean1),sigma=kronecker(temp,tXXi1))
          })),nrow=numpara)

          checksOC <- lapply(welk,function(h){
            Rorder <- as.matrix(RrO1[[h]][,-(1+numpara)])
            if(ncol(Rorder)==1){
              Rorder <- t(Rorder)
            }
            rorder <- as.matrix(RrO1[[h]][,1+numpara])
            apply(Rorder%*%randomDraws > rorder%*%t(rep(1,draws2)),2,prod)
          })
          relmeas <- rbind(relmeas,rep(1,2))
          relmeas[numhyp+1,] <- log(c(1,sum(Reduce("+",checksOC)==0)/draws2))
          rownames(relmeas)[numhyp+1] <- "complement"
        }
      }
    }
  }

  return(relmeas)
}

# The function computes the probability of an unconstrained draw falling in the complement subspace.
Student_prob_Hc <- function(mean1,scale1,df1,relmeas1,constraints,RrO1=NULL){

  numpara <- length(mean1)
  numhyp <- nrow(relmeas1)
  if(numhyp==1){
    relmeas <- t(relmeas1[1:numhyp,])
  }else{ relmeas <- relmeas1[1:numhyp,]}
  which_eq <- relmeas[,1] != log(1)
  if(sum(which_eq)==numhyp){ # Then the complement is equivalent to the unconstrained hypothesis.
    relmeas <- rbind(relmeas,log(rep(1,2)))
    rownames(relmeas)[numhyp+1] <- "complement"
  }else{ # So there is at least one hypothesis with only order constraints
    welk <- which(!which_eq)
    if(length(welk)==1){ # There is one hypothesis with only order constraints. Hc is complement of this hypothesis.
      relmeas <- rbind(relmeas,rep(1,2))
      relmeas[numhyp+1,2] <- 1 - exp(relmeas[welk,2])
      relmeas[numhyp+1,] <- log(relmeas[numhyp+1,])
      rownames(relmeas)[numhyp+1] <- "complement"
    }else{ # So more than one hypothesis with only order constraints
      # First we check whether ther is an overlap between the order constrained spaces.

      draws2 <- 1e4
      randomDraws <- rmvnorm(draws2,mean=rep(0,numpara),sigma=diag(numpara))
      #get draws that satisfy the constraints of the separate order constrained hypotheses
      checksOC <- lapply(welk,function(h){
        Rorder <- as.matrix(RrO1[[h]][,-(1+numpara)])
        if(ncol(Rorder)==1){
          Rorder <- t(Rorder)
        }
        rorder <- as.matrix(RrO1[[h]][,1+numpara])
        apply(randomDraws%*%t(Rorder) > rep(1,draws2)%*%t(rorder),1,prod)
      })
      checkOCplus <- Reduce("+",checksOC)

      if(sum(checkOCplus > 0) < draws2){ #then the joint order constrained hypotheses do not completely cover the parameter space.
        if(sum(checkOCplus>1)==0){ # then order constrained spaces are nonoverlapping
          relmeas <- rbind(relmeas,rep(log(1),2))
          relmeas[numhyp+1,2] <- 1 - sum(exp(relmeas[welk,2]))
          rownames(relmeas)[numhyp+1] <- "complement"
        }else{ #the order constrained subspaces at least partly overlap

          # the function below gives a rough estimate of the posterior probability under Hc
          # a bain type of algorithm but then using the Student distribution would be better.

          randomDraws <- rmvt(draws2,delta=mean1,sigma=scale1,df=df1)
          checksOC <- lapply(welk,function(h){
            Rorder <- as.matrix(RrO1[[h]][,-(1+numpara)])
            if(ncol(Rorder)==1){
              Rorder <- t(Rorder)
            }
            rorder <- as.matrix(RrO1[[h]][,1+numpara])
            apply(randomDraws%*%t(Rorder) > rep(1,draws2)%*%t(rorder),1,prod)
          })
          relmeas <- rbind(relmeas,rep(1,2))
          relmeas[numhyp+1,2] <- sum(Reduce("+",checksOC) == 0) / draws2
          relmeas[numhyp+1,] <- log(relmeas[numhyp+1,])
          rownames(relmeas)[numhyp+1] <- "complement"
        }
      }
    }
  }

  return(relmeas)
}

# from the output of the constraints in 'parse_hypothesis' create lists for the equality and order matrices
make_RrList <- function(parse_hyp){
  numhyp <- length(parse_hyp$hyp_mat)
  RrE <- lapply(1:numhyp,function(h){
    qE <- parse_hyp$n_constraints[h*2-1]
    if(qE==1){
      RrE_h <- t(as.matrix(parse_hyp$hyp_mat[[h]][1:qE,]))
    }else if(qE>1){
      RrE_h <- parse_hyp$hyp_mat[[h]][1:qE,]
    }else {RrE_h=NULL}
    RrE_h
  })
  RrO <- lapply(1:numhyp,function(h){
    qE <- parse_hyp$n_constraints[h*2-1]
    qO <- parse_hyp$n_constraints[h*2]
    if(qO==1){
      RrO_h <- t(as.matrix(parse_hyp$hyp_mat[[h]][qE+1:qO,]))
    }else if(qO>1){
      RrO_h <- parse_hyp$hyp_mat[[h]][qE+1:qO,]
    }else {RrO_h=NULL}
    RrO_h
  })
  return(list(RrE,RrO))
}

# from the output of the constraints in 'parse_hypothesis' create lists for the equality and order matrices
# different format parse_hyp object
make_RrList2 <- function(parse_hyp2){
  numhyp <- length(parse_hyp2$original_hypothesis)
  qE <- parse_hyp2$n_constraints[(0:(numhyp-1))*2+1]
  qO <- parse_hyp2$n_constraints[(1:numhyp)*2]
  RrE <- lapply(1:numhyp,function(h){
    startcon <- sum(qE[1:h]+qO[1:h])-qE[h]-qO[h]
    if(qE[h]==1){
      RrE_h <- t(as.matrix(parse_hyp2$hyp_mat[startcon+1:qE[h],]))
    }else if(qE[h]>1){
      RrE_h <- parse_hyp2$hyp_mat[startcon+1:qE[h],]
    }else {RrE_h=NULL}
    RrE_h
  })
  RrO <- lapply(1:numhyp,function(h){
    startcon <- sum(qE[1:h]+qO[1:h])-qE[h]-qO[h]
    if(qO[h]==1){
      RrO_h <- t(as.matrix(parse_hyp2$hyp_mat[startcon+qE[h]+1:qO[h],]))
    }else if(qO[h]>1){
      RrO_h <- parse_hyp2$hyp_mat[startcon+qE[h]+1:qO[h],]
    }else {RrO_h=NULL}
    RrO_h
  })
  return(list(RrE,RrO))
}

#for checking whether constraints are conflicting replace interval constraints by equality constraints
interval_RrStack <- function(RrStack){
  q1 <- nrow(RrStack)
  q2 <- ncol(RrStack)
  RrStack_out <- RrStack
  if(q1 > 1){
    row1 <- 1
    while(row1 < q1){
      for(row2 in (row1+1):q1){
        #        print(row2)
        if(sum(abs(RrStack_out[row1,-q2] + RrStack_out[row2,-q2]))==0){ # && RrStack_out[row1,q2]!=RrStack_out[row2,q2] ){
          #together row1 and row2 imply an interval constraint
          whichcol <- abs(RrStack_out[row1,-q2])!=0
          whichcol1 <- which(whichcol)
          if(sum(whichcol)==1){
            welkpos <- ifelse(RrStack_out[row1,c(whichcol,F)]>0,row1,row2)
            welkneg <- ifelse(RrStack_out[row1,c(whichcol,F)]<0,row1,row2)
            lb <- RrStack_out[welkpos,q2]
            ub <- -RrStack_out[welkneg,q2]
            RrStack_out[row1,] <- RrStack_out[welkpos,]
            RrStack_out[row1,q2] <- (ub+lb)/2
            RrStack_out <- RrStack_out[-row2,]
            q1 <- q1 - 1
          }else{
            RrStack_out[row1,q2] <- 0
            RrStack_out <- RrStack_out[-row2,]
            q1 <- q1 - 1
          }
          break
        }
      }
      row1 <- row1 + 1
    }
  }
  if(is.matrix(RrStack_out)==F){
    RrStack_out <- t(RrStack_out)
  }
  return(RrStack_out)
}

params_in_hyp <- function(hyp){
  params_in_hyp <- trimws(unique(strsplit(hyp, split = "[ =<>,\\(\\);&\\*+-]+", perl = TRUE)[[1]]))
  params_in_hyp <- params_in_hyp[!sapply(params_in_hyp, grepl, pattern = "^[0-9]*\\.?[0-9]+$")]
  params_in_hyp[grepl("^[a-zA-Z]", params_in_hyp)]
}


