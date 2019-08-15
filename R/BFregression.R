

#' @importFrom pracma rref
#' @importFrom mvtnorm dmvnorm pmvnorm rmvnorm dmvt pmvt rmvt
#' @importFrom Matrix rankMatrix
#' @importFrom stats rWishart qt
#' @importFrom MCMCpack rinvgamma
#' @importFrom MASS ginv
#' @method BF lm
#' @export
BF.lm <- function(x,
                  hypothesis = NULL,
                  prior = NULL,
                  parameter = NULL,
                  ...){

  # default BF on location parameters in a univariate normal linear model
  P <- 1
  N <- length(x$residuals)
  K <- length(x$coefficients) # dimension of predictors
  dummyX <- rep(F,K)
  names(dummyX) <- names(x$coefficients)
  Xmat <- model.matrix(x)
  Ymat <- model.matrix(x)%*%x$coefficients + x$residuals
  # check which are grouping covariates
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
    names_coef <- names(x$coefficients)
    # dummyX indicate which columns contain dummy group covariates
    dummyX1 <- apply(matrix(unlist(lapply(1:length(mains),function(faclev){
      unlist(lapply(1:length(names_coef),function(cf){
        grepl(mains[faclev],names_coef[cf])
      }))
    })),nrow=length(names_coef)),1,max)==1
    if(is.matrix(apply(Xmat,2,table))){
      if(nrow(apply(Xmat,2,table))){
        dummyX2 <- rep(T,ncol(Xmat))
      }else{
        dummyX2 <- rep(F,ncol(Xmat))
      }
    }else{
      dummyX2 <- unlist(lapply(apply(Xmat,2,table),length))==2
    }
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
  names_coef <- names(x$coefficients)

  # prior hyperparameters
  df0 <- 1 # should be the same as sum(rep(bj,times=Nj))-K-P+1
  Scale0 <- kronecker(S_b,tXXi_b)
  mean0 <- as.matrix(rep(0,K*P))
  # posterior hyperparameters
  dfN <- N-K-P+1
  ScaleN <- kronecker(S,tXXi)/(N-K-P+1) # off-diagonal elements have no meaning
  meanN <- as.matrix(c(BetaHat))

  # Hypotheses for exploratory test
  # H0: beta = 0
  # H1: beta < 0
  # H2: beta > 0
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

  BFtu_exploratory <- relfit / relcomp
  colnames(BFtu_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
  PHP_exploratory <- BFtu_exploratory / apply(BFtu_exploratory,1,sum)

  # Additional exploratory tests in the case of an aov-object
  if(class(x)[1]=="aov"){
    # check main effects
    BFmain <- unlist(lapply(1:length(numlevels),function(fac){
      name1 <- names(numlevels[fac])
      mains1 <- mains[sum(numlevels[1:fac])-numlevels[fac]+1:numlevels[fac]]
      which0 <- unlist(lapply(1:length(colnames(Xmat)),function(col){
        sum(colnames(Xmat)[col]==mains1)==1
      }))
      if(sum(which0)>0){
        RrE_f <- matrix(0,nrow=sum(which0),ncol=length(colnames(Xmat))+1)
        for(r1 in 1:sum(which0)){RrE_f[r1,which(which0)[r1]]<-1}
        relcomp_f <- Student_measures(mean1=mean0,Scale1=Scale0,df1=df0,RrE1=RrE_f,RrO1=NULL)
        relfit_f <- Student_measures(mean1=meanN,Scale1=ScaleN,df1=dfN,RrE1=RrE_f,RrO1=NULL)
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
    for(c1 in 1:ncol(matcov)){
      if(c1 < ncol(matcov)){
        numeffects_c <- sum(matcov[,c1])
        if(numeffects_c>1){
          count_interaction <- count_interaction + 1
          interactionset <- names(which(matcov[,c1]))
          whichx <- apply(matcov[which(matcov[,c1]),],2,sum)==length(interactionset)
          RrE_ia <- matrix(0,nrow=sum(whichx),ncol=K+1)
          for(r1 in 1:sum(whichx)){RrE_ia[r1,which(whichx)[r1]]<-1}
          relcomp_ia <- Student_measures(mean1=mean0,Scale1=Scale0,df1=df0,RrE1=RrE_ia,RrO1=NULL)
          names(relcomp_ia) <- c("c=","c>")
          relfit_ia <- Student_measures(mean1=meanN,Scale1=ScaleN,df1=dfN,RrE1=RrE_ia,RrO1=NULL)
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

  # compute posterior estimates
  postestimates <- cbind(meanN,meanN,
                         t(matrix(unlist(lapply(1:length(meanN),function(coef){
                           ub <- qt(p=.975,df=dfN)*sqrt(ScaleN[coef,coef])+meanN[coef,1]
                           lb <- qt(p=.025,df=dfN)*sqrt(ScaleN[coef,coef])+meanN[coef,1]
                           return(c(lb,ub))
                         })),nrow=2))
  )
  row.names(postestimates) <- names_coef
  colnames(postestimates) <- c("mean","median","2.5%","97.5%")

  # confirmatory BF test
  if(!is.null(hypothesis)){
    #read constraints
    parse_hyp <- parse_hypothesis(names_coef,hypothesis)
    RrList <- make_RrList2(parse_hyp)
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]

    # RrStack is used to check conflicting constraints, and for the default prior location
    RrStack <- rbind(do.call(rbind,RrE),do.call(rbind,RrO))
    RrStack <- interval_RrStack(RrStack)
    if(nrow(RrStack)>1){
      RStack <- RrStack[,-(K+1)]
      rStack <- RrStack[,(K+1)]
    }else{
      RStack <- matrix(RrStack[,-(K+1)],nrow=1)
      rStack <- RrStack[,(K+1)]
    }

    # check if a common boundary exists for prior location under all constrained hypotheses
    if(nrow(RrStack) > 1){
      rref_ei <- rref(RrStack)
      nonzero <- RrStack[,K+1]!=0
      if(max(nonzero)>0){
        row1 <- max(which(nonzero))
        if(sum(abs(rref_ei[row1,1:K]))==0){
          stop("No common boundary point for prior location. Conflicting constraints.")
        }
      }
    }
    # default prior location
    adjmean <- ginv(RStack)%*%rStack

    # prior hyperparameters
    df0 <- 1 # should be the same as sum(rep(bj,times=Nj))-K-P+1
    Scale0 <- kronecker(S_b,tXXi_b)
    mean0 <- adjmean
    # posterior hyperparameters
    dfN <- N-K-P+1
    ScaleN <- kronecker(S,tXXi)/(N-K-P+1) # off-diagonal elements have no meaning
    meanN <- as.matrix(c(BetaHat))

    #number of hypotheses that are specified
    numhyp <- length(RrO)

    relcomp <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Student_measures(mean1=mean0,Scale1=Scale0,df1=df0,RrE1=RrE[[h]],RrO1=RrO[[h]],
                       names1=names_coef,constraints1=parse_hyp$original_hypothesis[h])
    })),nrow=2))

    relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Student_measures(meanN,ScaleN,dfN,RrE[[h]],RrO[[h]],
                       names1=names_coef,constraints1=parse_hyp$original_hypothesis[h])
    })),nrow=2))

    # Compute relative fit/complexity for the complement hypothesis
    relfit <- Student_prob_Hc(meanN,ScaleN,dfN,relfit,hypothesis,RrO)
    relcomp <- Student_prob_Hc(mean0,Scale0,df0,relcomp,hypothesis,RrO)
    row.names(relcomp)[1:numhyp] <- parse_hyp$original_hypothesis
    row.names(relfit)[1:numhyp] <- parse_hyp$original_hypothesis
    colnames(relcomp) <- c("c_E","c_O")
    colnames(relfit) <- c("f_E","f_O")

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
    #tested hypotheses
    hypotheses <- row.names(relfit)
  }else{
    BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relfit <-
      relcomp <- hypotheses <- BFtable <- priorprobs <- NULL
  }

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
    call=match.call())

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
  relE <- relO <- 1
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
    relE <- mean(unlist(lapply(covm1_E,function(temp) dmvnorm(mean1_E,mean=mean1_E,sigma=temp))))

  }else{
    if(is.null(RrE1) && !is.null(RrO1)){ #only order constraints
      RO1 <- RrO1[,-(K*P+1)]
      if(!is.matrix(RO1)){
        RO1 <- matrix(RO1,ncol=K*P)
      }
      qO1 <- nrow(RO1)
      rO1 <- RrO1[,(K*P+1)]

      if(rankMatrix(RO1)[[1]]==nrow(RO1)){ #RO1 is of full row rank. So use transformation.

        Scale1inv <- solve(Scale1)
        relO <- mean(unlist(lapply(1:1e3,function(s){
          Sigma1 <- solve(rWishart(1,df=df1+P-1,Sigma=Scale1inv)[,,1])
          meanO <- c(RO1%*%mean1)
          covmO <- RO1%*%kronecker(Sigma1,tXXi1)%*%t(RO1)
          pmvnorm(lower=rO1,upper=Inf,mean=meanO,sigma=covmO)
        })))

      }else{ #no linear transformation can be used; pmvt cannot be used. Use bain with a multivariate normal approximation
        #compute covariance matrix for multivariate normal distribution
        mean1 <- c(Mean1)
        names(mean1) <- c(Names1)
        if(df1>2){ #posterior measures
          covm1 <- kronecker(Scale1,tXXi1)/(df1-2)
          bain_res <- bain(x=mean1,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
          relO <- bain_res$fit[1,3]
        }else if(df1==2){ #posterior measures
          covm1 <- kronecker(Scale1,tXXi1)/(df1-1)
          bain_res <- bain(x=mean1,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
          relO <- bain_res$fit[1,3]
        }else{
          covm1 <- kronecker(Scale1,tXXi1) #for prior with df1==1, probability independent of common factor of scale1
          bain_res <- bain(x=mean1,hypothesis=constraints1,Sigma=covm1,n=df1) #n not used in computation
          relO <- bain_res$fit[1,4]
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
      relE <- mean(unlist(lapply(covm1_E,function(temp) dmvnorm(c(mean1_E),mean=mean1_E,sigma=temp))))

      if(rankMatrix(RrO1)[[1]] == nrow(RrO1)){
        covm1_O <- lapply(SigmaList,function(temp) R1%*%(kronecker(temp,tXXi1))%*%t(R1) )
        mean1_O <- c(RO1%*%mean1) - rO1  ## mu_zeta_O in the paper

        mean1_OE <- lapply(covm1_O,function(temp){
          as.vector(mean1_O +
                      matrix(temp[(qE1+1):qC1,1:qE1],ncol=qE1)%*%solve(temp[1:qE1,1:qE1])%*%(rE1-mean1_E))
        })
        covm1_OE <- lapply(covm1_O,function(temp) temp[(qE1+1):qC1,(qE1+1):qC1] -
                             matrix(temp[(qE1+1):qC1,1:qE1],ncol=qE1)%*%solve(temp[1:qE1,1:qE1])%*%
                             matrix(temp[1:qE1,(qE1+1):qC1],nrow=qE1))
        #check covariance because some can be nonsymmetric due to a generation error
        welk1 <- which(unlist(lapply(covm1_OE,function(temp) isSymmetric(temp,
                             tol = sqrt(.Machine$double.eps),check.attributes = FALSE)))==T)
        covm1_OE <- covm1_OE[welk1]
        mean1_OE <- mean1_OE[welk1]
        relO <- mean(mapply(function(mu_temp,Sigma_temp) pmvnorm(lower=rep(0,qO1),
                             upper=rep(Inf,qO1),mean=mu_temp,sigma=Sigma_temp)[1],mean1_OE,covm1_OE))
      }else{ #use bain for the computation of the probability

        mean1 <- c(Mean1)
        names(mean1) <- c(Names1)
        if(df1>2){ #posterior measures
          covm1 <- kronecker(Scale1,tXXi1)/(df1-2)
          bain_res <- bain(x=mean1,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
          relO <- bain_res$fit[1,3]
        }else if(df1==2){ #posterior measures
          covm1 <- kronecker(Scale1,tXXi1)/(df1-1)
          bain_res <- bain(x=mean1,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
          relO <- bain_res$fit[1,3]
        }else{
          covm1 <- kronecker(Scale1,tXXi1) #for prior with df1==1, probability independent of common factor of scale1
          bain_res <- bain(x=mean1,hypothesis=constraints1,Sigma=covm1,n=df1) #n not used in computation
          relO <- bain_res$fit[1,4]
        }
      }
    }
  }
  return(c(relE,relO))
}

# compute relative meausures (fit or complexity) under a multivariate Student t distribution
Student_measures <- function(mean1,Scale1,df1,RrE1,RrO1,names1=NULL,constraints1=NULL){ # Volgens mij moet je hier ook N meegeven
  K <- length(mean1)
  relE <- relO <- 1
  if(!is.null(RrE1) && is.null(RrO1)){ #only equality constraints
    RE1 <- RrE1[,-(K+1)]
    if(!is.matrix(RE1)){
      RE1 <- matrix(RE1,ncol=K)
    }
    rE1 <- RrE1[,(K+1)]
    qE1 <- nrow(RE1)
    meanE <- RE1%*%mean1
    scaleE <- RE1%*%Scale1%*%t(RE1)
    relE <- dmvt(rE1,delta=c(meanE),sigma=scaleE,df=df1,log=FALSE)
  }
  if(is.null(RrE1) && !is.null(RrO1)){ #only order constraints
    RO1 <- RrO1[,-(K+1)]
    if(!is.matrix(RO1)){
      RO1 <- matrix(RO1,ncol=K)
    }
    qO1 <- nrow(RO1)
    rO1 <- RrO1[,(K+1)]

    if(rankMatrix(RO1)[[1]]==nrow(RO1)){ #RO1 is of full row rank. So use transformation.
      meanO <- c(RO1%*%mean1)
      scaleO <- RO1%*%Scale1%*%t(RO1)
      relO <- ifelse(nrow(scaleO)==1,
                     pt((rO1-meanO)/sqrt(scaleO[1,1]),df=df1,lower.tail=FALSE), #univariate
                     pmvt(lower=rO1,upper=Inf,delta=meanO,sigma=scaleO,df=df1,type="shifted")) #multivariate
    }else{ #no linear transformation can be used; pmvt cannot be used. Use bain with a multivariate normal approximation
      #compute covariance matrix for multivariate normal distribution
      row.names(mean1) <- names1
      if(df1>2){ # we need posterior measures
        covm1 <- Scale1/(df1-2)
        mean1vec <- c(mean1)
        names(mean1vec) <- row.names(mean1)
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
        relO <- bain_res$fit[1,3]
      }else if(df1==2){ # we need posterior measures (there is very little information)
        covm1 <- Scale1/(df1-1)
        mean1vec <- c(mean1)
        names(mean1vec) <- row.names(mean1)
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
        relO <- bain_res$fit[1,3]
      }else{ #then df=1, so we need prior measures
        covm1 <- Scale1 #for prior with df1==1, probability independent of common factor of scale1
        mean1vec <- c(mean1)
        names(mean1vec) <- row.names(mean1)
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=covm1,n=df1) #n not used in computation
        relO <- bain_res$fit[1,4]
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
    Tm <- rbind(RE1, D2)

    #b)
    Tmean1 <- Tm %*% mean1
    Tscale1 <- Tm %*% Scale1 %*% t(Tm)

    # relative meausure for equalities
    relE <- dmvt(x = t(rE1), delta = Tmean1[1:qE1], sigma = matrix(Tscale1[1:qE1, 1:qE1], ncol = qE1), df = df1, log = FALSE)

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

    if(rankMatrix(RO1tilde)[[1]] == nrow(RO1tilde)){
      rO1tilde <- as.vector(rO1tilde)

      delta_trans <- as.vector(RO1tilde %*% Tmean1OgE)
      scale1_trans <- RO1tilde %*% Tscale1OgE %*% t(RO1tilde)

      if(nrow(scale1_trans) == 1){ # univariate
        relO <- pt((rO1tilde - delta_trans) / sqrt(scale1_trans), df = df1+qE1, lower.tail = FALSE)[1]
      } else { # multivariate
        relO <- pmvt(lower = rO1tilde, upper = Inf, delta = delta_trans, sigma = scale1_trans, df = df1+qE1, type = "shifted")[1]
      }

    }else{ #use bain for the computation of the probability
      #compute covariance matrix for multivariate normal distribution
      row.names(mean1) <- names1
      if(df1>2){ # we need posterior measures
        covm1 <- Scale1/(df1-2)
        mean1vec <- c(mean1)
        names(mean1vec) <- row.names(mean1)
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
        relO <- bain_res$fit[1,3]
      }else if(df1==2){ # we need posterior measures (there is very little information)
        covm1 <- Scale1/(df1-1)
        mean1vec <- c(mean1)
        names(mean1vec) <- row.names(mean1)
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=covm1,n=999) #n not used in computation
        relO <- bain_res$fit[1,3]
      }else{ #then df=1, so we need prior measures
        covm1 <- Scale1 #for prior with df1==1, probability independent of common factor of scale1
        mean1vec <- c(mean1)
        names(mean1vec) <- row.names(mean1)
        bain_res <- bain(x=mean1vec,hypothesis=constraints1,Sigma=covm1,n=df1) #n not used in computation
        relO <- bain_res$fit[1,4]
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
  which_eq <- relmeas[,1] != 1
  if(sum(which_eq)==numhyp){ # Then the complement is equivalent to the unconstrained hypothesis.
    relmeas <- rbind(relmeas,rep(1,2))
    rownames(relmeas)[numhyp+1] <- "complement"
  }else{ # So there is at least one hypothesis with only order constraints
    welk <- which(!which_eq)
    if(length(welk)==1){ # There is one hypothesis with only order constraints. Hc is complement of this hypothesis.
      relmeas <- rbind(relmeas,rep(1,2))
      relmeas[numhyp+1,2] <- 1 - relmeas[welk,2]
      rownames(relmeas)[numhyp+1] <- "complement"
    }else{ # So more than one hypothesis with only order constraints
      # First we check whether ther is an overlap between the order constrained spaces.

      # Caspar, here we need the RE and RO which are lists of
      # matrices for equality and order constraints under the hypotheses. We can probably do this
      # using the R-code you wrote and a vector of names of the correlations but I don't know
      # how exactly. When running your function I also get an error message saying that he
      # does not know the function "rename_function".

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
          relmeas[numhyp+1,2] <- 1 - sum(relmeas[welk,2])
          rownames(relmeas)[numhyp+1] <- "complement"
        }else{ #the order constrained subspaces at least partly overlap

          # funtion below gives a rough estimate of the posterior probability under Hc
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
          relmeas[numhyp+1,] <- c(1,sum(Reduce("+",checksOC)==0)/draws2)
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
  which_eq <- relmeas[,1] != 1
  if(sum(which_eq)==numhyp){ # Then the complement is equivalent to the unconstrained hypothesis.
    relmeas <- rbind(relmeas,rep(1,2))
    rownames(relmeas)[numhyp+1] <- "complement"
  }else{ # So there is at least one hypothesis with only order constraints
    welk <- which(!which_eq)
    if(length(welk)==1){ # There is one hypothesis with only order constraints. Hc is complement of this hypothesis.
      relmeas <- rbind(relmeas,rep(1,2))
      relmeas[numhyp+1,2] <- 1 - relmeas[welk,2]
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
          relmeas[numhyp+1,2] <- 1 - sum(relmeas[welk,2])
          rownames(relmeas)[numhyp+1] <- "complement"
        }else{ #the order constrained subspaces at least partly overlap

          # funtion below gives a rough estimate of the posterior probability under Hc
          # a bain type of algorithm would be better of course. but for now this is ok.

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
# RrStack <- rbind(do.call(rbind,RrE),do.call(rbind,RrO))
# row.names(RrStack) <- colnames(RrStack) <- NULL
interval_RrStack <- function(RrStack){
  q1 <- nrow(RrStack)
  q2 <- ncol(RrStack)
  RrStack_out <- RrStack
  if(q1 > 1){
    row1 <- 1
    while(row1 < q1){
      for(row2 in (row1+1):q1){
#        print(row2)
        if(sum(abs(RrStack_out[row1,-q2] + RrStack_out[row2,-q2]))==0 && RrStack_out[row1,q2]!=RrStack_out[row2,q2] ){
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





