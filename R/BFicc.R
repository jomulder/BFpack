### Joris Mulder 2019. Bayes factor testing of multiple random intercept models
### via multiple lmer-objects based on Mulder & Fox (2013, 2019).
### with extension to unbalanced data

#' @importFrom stats rgamma rnorm rbeta dbeta runif
#' @importFrom lme4 getME VarCorr
#' @method BF lmerMod
#' @export
BF.lmerMod <- function(x,
                           hypothesis = NULL,
                           prior = NULL,
                           complement = TRUE,
                           ...){

  numcat <- length(x@cnms)
  namescat <- unlist(lapply(1:numcat,function(ca){
    x@cnms[[ca]]
  }))
  if(numcat>1){
    iccnames <- unlist(lapply(1:numcat,function(nc){namescat[nc]}))
  }else{ iccnames <- "icc" }

  # check if the lmer-model only has a random intercept or category specific random intercepts
  Zstack <- Reduce(cbind,getME(x,"mmList"))
  if(numcat>1){ #check if the random effects are category specific random intercepts
    for(ca in 1:numcat){
      freq_ca <- table(Zstack[,ca])
      if(sum(abs(sort(as.integer(names(freq_ca))) - c(0,1))) !=0){
        stop("only models with a single random intercept or category specific random intercepts are currently supported when testing icc's.")
      }
    }
  }else{
    freq_ca <- table(getME(x,"mmList")[[1]])
    if(as.integer(names(freq_ca))!=1 || length(as.integer(names(freq_ca)))!=1){
      stop("only models with a single random intercept or category specific random intercepts are currently supported when testing icc's.")
    }
  }

  # sort data per cluster
  clusterindex <- x@flist[[1]]
  if(length(table(table(clusterindex)))>1){
    balanced <- FALSE
    pvec <- table(clusterindex)
  }else{
    balanced <- TRUE
    pvec <- rep(table(clusterindex)[1],length(table(clusterindex)))
  }

  levels(clusterindex) <- 1:length(levels(clusterindex))
  reorder1 <- order(as.integer(clusterindex))
  nclusters <- length(levels(clusterindex)) #total number of groups/clusters
  ystack <- getME(x,"y")[reorder1]
  Xstack <- as.matrix(getME(x,"X")[reorder1,])
  Zstack <- Zstack[reorder1,]
  # next sort data per category
  if(numcat>1){
    # sort rows per category
    firstofcluster <- cumsum(c(0,pvec[1:(nclusters-1)])) + 1
    catassign <- unlist(lapply(firstofcluster,function(cluster){
      ifelse(sum(Zstack[cluster,1:numcat])==1,
             which(Zstack[cluster,1:numcat]==1),
             0)
    }))
    if(sum(names(table(catassign))=="0")==0){ #all groups belong to a category
      ngroups <- table(catassign)
      numgroups <- sum(ngroups)
      names(ngroups) <- NULL
      #reorder data matrices according to categories
      reorder2 <- unlist(lapply(1:numcat,function(ca){
        welk_ca <- which(catassign==ca)
        names(welk_ca) <- NULL
        unlist(lapply(welk_ca,function(cluster){
          firstofcluster[cluster]+0:(pvec[cluster]-1)
        }))
      }))
      #update order of group sizes
      pvec <- pvec[unlist(lapply(1:numcat,function(ca){
        welk_ca <- which(catassign==ca)
      }))]
      Xstack <- Xstack[reorder2,]
      ystack <- ystack[reorder2]
      Zstack <- Zstack[reorder2,]
    }else{ #only include groups that belong to a category
      stop("Some groups don't belong to a group category. Exclude these groups from the data.")
    }
  }else{
    ngroups <- numgroups <- nclusters
  }

  #transform data matrices with Helmert matrix
  firstofcluster0 <- c(0,cumsum(pvec[1:(numgroups-1)]))
  zWstack <- do.call(rbind,
                     lapply(1:length(pvec),function(j){
                       H_j <- Helmert(pvec[j])
                       H_j%*%cbind(ystack[firstofcluster0[j]+1:pvec[j]],
                                   as.matrix(Xstack[firstofcluster0[j]+1:pvec[j],]) )
                     })
  )
  #extract ML estimates for rho
  tau2ML <- unlist(VarCorr(x))
  sigma2ML <- attr(VarCorr(x),"sc")**2
  rhoML <- tau2ML/(tau2ML+sigma2ML)

  shape0 <- c(1,1) # set uniform priors for icc's
  cat("First, unconstrained icc analysis...")
  cat("\n")
  cat("\n")
  numuncdraws <- 5e4
  marglike_Hu <- MargLikeICC_Hq(rhoML,zW=zWstack,ngroups,pvec,samsize1=numuncdraws,
                                samsize2=4e4,unique1=1:numcat)

  postestimates <- marglike_Hu[[4]]
  colnames(postestimates) <- iccnames
  cat("Second, exploratory testing of icc's...")
  cat("\n")
  BFtu_exploratory_icc <- t(matrix(unlist(lapply(1:numcat,function(nc){

    cat(paste0(iccnames[nc],"; "))

    marglike_explo <- rep(0,3)

    if(numcat>1){
      unique_c <- rep(1,numcat)
      unique_c[nc] <- 0
      unique_c[-nc] <- 1:(numcat-1)
    }else {
      unique_c <- 0
    }

    # zero icc
    marglike_explo[1] <- MargLikeICC_Hq(rhoML,zWstack,ngroups,pvec,unique1=unique_c)[[1]]

    # positive icc
    marglike_positive <- marglike_Hu[[1]] + log(marglike_Hu$postprobpositive[nc]) -
      log(marglike_Hu$priorprobpositive[nc])
    marglike_explo[3] <- marglike_positive
    # negative icc
    marglike_negative <- marglike_Hu[[1]] + log(1-marglike_Hu$postprobpositive[nc]) -
      log(1-marglike_Hu$priorprobpositive[nc])
    marglike_explo[2] <- marglike_negative

    return(exp(marglike_explo - marglike_Hu[[1]]))
  })),nrow=3))

  colnames(BFtu_exploratory_icc) <- c("icc=0","icc<0","icc>0")
  row.names(BFtu_exploratory_icc) <- iccnames
  PHP_exploratory_icc <- round(BFtu_exploratory_icc / apply(BFtu_exploratory_icc,1,sum),3)
  priorprobs <- rep(1,3)/3 #prior probs for exploratory tests
  cat("\n")
  cat("\n")
  if(!is.null(hypothesis)){ # confirmatory test with constrained hypotheses on icc's.
    cat("Third, confirmatory testing of icc's...")
    cat("\n")
    parse_hyp <- parse_hypothesis(iccnames,hypothesis)
    parse_hyp$hyp_mat <- do.call(rbind, parse_hyp$hyp_mat)
    RrList <- make_RrList2(parse_hyp)
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]
    # check if icc's are only tested against each other or against zero
    numhyp <- length(RrE)
    for(h in 1:numhyp){
      if(!is.null(RrE[[h]])){
        for(r in 1:nrow(RrE[[h]])){
          row1 <- RrE[[h]][r,]
          if( !(sum(abs(row1))==1 || sum(row1)==0) ){
            stop("icc's can only be compared with each other or to zero.")
          }
        }
      }
      if(!is.null(RrO[[h]])){
        for(r in 1:nrow(RrO[[h]])){
          freq1 <- table(sort(RrO[[h]][r,]))
          row1 <- RrO[[h]][r,]
          if( !(sum(abs(row1))==1 || sum(row1)==0) ){
            stop("icc's can only be compared with each other or to zero.")
          }
        }
      }
    }

    output_marglike_icc <- t(matrix(unlist(lapply(1:numhyp, function(h){

      cat(paste0(parse_hyp$original_hypothesis[h],"; "))
      cat("\n")

      # get unconstrained prior draws, if needed for computing prior probabilities
      # in case of only order constraints
      pcat <- rep(1:length(ngroups),times=ngroups)
      unique1 <- 1:numcat
      LB <- unlist(unlist(lapply(1:numcat,function(c){
        -1/(max(pvec[pcat==c])-1)
      })))
      priordraws <- matrix(unlist(lapply(1:numcat,function(c){
        rbeta(numuncdraws,1,1) * (1 - LB[c]) + LB[c]
      })),ncol=numcat)

      # code equal icc's with same integer for marglike2_Hq function
      unique_h <- 1:numcat
      if(!is.null(RrE[[h]])){
        unique_h <- rep(NA,numcat)
        zeroMat <- RrE[[h]][which(apply(RrE[[h]],1,sum)==1),]
        if(!is.matrix(zeroMat)){
          zeroMat <- matrix(zeroMat,nrow=1)
        }
        unique_h[which(apply(zeroMat,2,sum)==1)] <- 0
        teller <- 0
        for(row1 in which(apply(RrE[[h]],1,sum)==0)){
          welk1 <- which(RrE[[h]][row1,]!=0)
          isna_h <- is.na(unique_h[welk1])
          if(sum(isna_h)==2){
            teller <- teller + 1
            unique_h[welk1] <- teller
          }else{ #one is already assigned a unique code
            unique_h[welk1] <- unique_h[welk1[!isna_h]]
            #
          }

        }
        if(sum(is.na(unique_h))>0){ #unconstrained icc's receive unique code
          unique_h[is.na(unique_h)] <- teller + 1:sum(is.na(unique_h))
          teller <- teller + sum(is.na(unique_h))
        }
      }
      if(!is.null(RrO[[h]])){
        unicum <- unique(unique_h[unique_h!=0])
        inequalities_h <- matrix(0,nrow=nrow(RrO[[h]]),ncol=max(unicum)+1)
        for(u in sort(unicum)){

          welk <- which(unique_h == u)
          if(length(welk) > 1){
            if(nrow(RrO[[h]]) > 1){
              inequalities_h[,u] <- apply(RrO[[h]][,which(unique_h == u)],1,sum)
            }else{
              inequalities_h[,u] <- apply(t(RrO[[h]][,which(unique_h == u)]),1,sum)
            }
          }else{ #length is 1
            inequalities_h[,u] <- RrO[[h]][,which(unique_h == u)]
          }
        }
      } else inequalities_h = 0
      if(is.null(RrE[[h]])){ #only order constraints; use output from unconstrained analysis
        priorprob_h <- mean(apply(cbind(priordraws,rep(-1,length(numuncdraws)))%*%t(inequalities_h)>0,1,prod))
        postprob_h <- mean(apply(cbind(marglike_Hu$postdraws,rep(-1,length(numuncdraws)))%*%
                                   t(inequalities_h)>0,1,prod))
        marglike2_h <- list(marglike_Hu[[1]],postprob_h,priorprob_h)
      }else{
        marglike2_h <- MargLikeICC_Hq(rhoML,zW=zWstack,ngroups,pvec,unique1=unique_h,
                                      inequalities=inequalities_h)[1:3]
      }
      return(c(unlist(marglike2_h),ifelse(is.null(RrE[[h]]),1,0)))
    })),nrow=4))
    if(complement == TRUE){
      #compute BF for complement hypothesis
      if(sum(output_marglike_icc[,4])==0){ #the complement is equivalent to the unconstrained model
        output_marglike_icc <- rbind(output_marglike_icc,c(unlist(marglike_Hu)[1:3],1))
      } else { #the complement is the complement of the joint of the order hypotheses
        which_order <- which(output_marglike_icc[,4]==1)
        if(length(which_order)==1){
          probs <- output_marglike_icc[which_order,2:3]
          marglike_Hc <- marglike_Hu[[1]] + log(1-probs[1]) - log(1-probs[2])
          output_marglike_icc <- rbind(output_marglike_icc,c(marglike_Hc,1-probs[1],1-probs[2],1))
        }else {
          probs <- apply(output_marglike_icc[which_order,2:3],2,sum)
          marglike_Hc <- marglike_Hu[[1]] + log(1-probs[1]) - log(1-probs[2])
          output_marglike_icc <- rbind(output_marglike_icc,c(marglike_Hc,1-probs[1],1-probs[2],1))
        }
      }
      row.names(output_marglike_icc) <- c(parse_hyp$original_hypothesis,"complement")
    }else{
      row.names(output_marglike_icc) <- c(parse_hyp$original_hypothesis)
    }
    relcomp <- matrix(c(rep(NA,nrow(output_marglike_icc)),output_marglike_icc[,3]),ncol=2)
    relfit <- matrix(c(rep(NA,nrow(output_marglike_icc)),output_marglike_icc[,2]),ncol=2)
    #compute log marginal likelihood for H* without order constraints
    BF_E <- exp(output_marglike_icc[,1] - log(output_marglike_icc[,2]) + log(output_marglike_icc[,3]) - marglike_Hu[[1]])
    BFtu_confirmatory_icc <- exp(output_marglike_icc[,1] - marglike_Hu[[1]])
    #compute BFmatrix and PHPs
    logBFmatrix <- matrix(rep(output_marglike_icc[,1],numhyp+complement),nrow=numhyp+complement) -
      matrix(rep(output_marglike_icc[,1],each=numhyp+complement),nrow=numhyp+complement)
    if(complement == TRUE){
      row.names(logBFmatrix) <- colnames(logBFmatrix) <- c(parse_hyp$original_hypothesis,"complement")
    }else{
      row.names(logBFmatrix) <- colnames(logBFmatrix) <- c(parse_hyp$original_hypothesis)
    }
    BFmatrix_confirmatory_icc <- round(exp(logBFmatrix),3)
    BFta_confirmatory_icc <- exp(output_marglike_icc[,1] - max(output_marglike_icc[,1]))
    # Change prior probs in case of default setting
    if(is.null(prior)){
      priorprobs <- rep(1/length(BFtu_confirmatory_icc),length(BFtu_confirmatory_icc))
    }else{
      if(!is.numeric(prior) || length(prior)!=length(BFtu_confirmatory_icc)){
        warning(paste0("Argument 'prior' should be numeric and of length ",as.character(length(BFtu_confirmatory_icc)),". Equal prior probabilities are used."))
        priorprobs <- rep(1/length(BFtu_confirmatory_icc),length(BFtu_confirmatory_icc))
      }else{
        priorprobs <- prior
      }
    }
    PHP_confirmatory_icc <- priorprobs*BFta_confirmatory_icc / sum(priorprobs*BFta_confirmatory_icc)
    BFtable <- cbind(relcomp,relfit,BF_E,relfit[,2]/relcomp[,2],
                     BF_E*relfit[,2]/relcomp[,2],PHP_confirmatory_icc)
    row.names(BFtable) <- names(PHP_confirmatory_icc)
    colnames(BFtable) <- c("comp_E","comp_O","fit_E","fit_O","BF_E","BF_O","BF","PHP")
    hypotheses <- names(BFta_confirmatory_icc)

  }else{
    BFmatrix_confirmatory_icc <- PHP_confirmatory_icc <- BFtu_confirmatory_icc <- relfit <-
      relcomp <- hypotheses <- BFtable <- priorprobs <- NULL
  }
  BFlm_out <- list(
    BFtu_exploratory=BFtu_exploratory_icc,
    PHP_exploratory=PHP_exploratory_icc,
    BFtu_confirmatory=BFtu_confirmatory_icc,
    PHP_confirmatory=PHP_confirmatory_icc,
    BFmatrix_confirmatory=BFmatrix_confirmatory_icc,
    BFtable_confirmatory=BFtable,
    prior=priorprobs,
    hypotheses=hypotheses,
    estimates=postestimates,
    model=x,
    bayesfactor="Bayes factor based on uniform priors",
    parameter="intraclass correlations",
    call=match.call())

  class(BFlm_out) <- "BF"

  return(BFlm_out)

}

int_lhood <- function(rhoS,ngroups,pvec,N,K,Wmat1,zvec1,tWW2,tWz2,tzz2){
  #the log integrated likelihood

  diagDi <- c( (1-rep(rhoS,times=ngroups))/(1+(pvec-1)*rep(rhoS,times=ngroups)) )
  tWDiW1 <- t(Wmat1)%*%(diagDi*Wmat1)
  tWDiz1 <- t(Wmat1)%*%(diagDi*zvec1)
  tzDiz1 <- sum(zvec1**2*diagDi)
  s2 <- tzDiz1+tzz2 - c(t(tWDiz1+tWz2)%*%solve(tWDiW1+tWW2)%*%(tWDiz1+tWz2) )
  return(
    .5*sum(log(diagDi)) - .5*log(det(tWDiW1+tWW2)) - (N-K)/2*log(s2)
  )
}

#' @importFrom stats rnorm rbeta dbeta
MargLikeICC_Hq <- function(rhoML,zW,ngroups,pvec,samsize1=4e4,samsize2=4e4,
                           unique1,inequalities=0,complement=FALSE){

  #E.g., for categories=5, unique1=c(0,2,0,1,1),inequalities=[1 -1 0],
  #the hypothesis equals, Hq:rho1=rho3=0, rho4=rho5>rho2
  #samsize1 sets the sample size for the number of draws from the proposal distribution for the importance
  #sample estimate of the marginal likelihood excluding the inequality constraints.
  #samsize2 sets the sample size for computing the probability that the inequality constraints hold,
  #and for sampling from the proposal distribution. to get the IS estimate.

  numcat <- length(ngroups)
  N <- sum(pvec)

  zvec <- zW[,1]
  Wmat <- as.matrix(zW[,-1])

  K <- ncol(Wmat)
  numgroups <- sum(ngroups)
  firstofcluster <- c(0,cumsum(pvec[1:(numgroups-1)]))+1
  restcluster <- (1:N)[-firstofcluster]
  zvec1 <- zvec[firstofcluster]
  Wmat1 <- as.matrix(Wmat[firstofcluster,])
  zvec2 <- zvec[restcluster]
  Wmat2 <- as.matrix(Wmat[restcluster,])
  tWW2 <- t(Wmat2)%*%Wmat2
  tWz2 <- t(Wmat2)%*%zvec2
  tzz2 <- sum(zvec2**2)

  if(sum(unique1)==0){ #null model with no random effects
    marglike <- int_lhood(rhoS=rep(0,length=numcat),ngroups,pvec,N,K,Wmat1,zvec1,tWW2,tWz2,tzz2) +
      lgamma((N-K)/2) - (N-K)/2*log(pi)
    postprobpositive <- priorprobpositive <- postestimates <- postdraws <- drawsMat <- NULL
    priorprob <- postprob <- 1

  }else{
    numHrho <- max(unique1)

    transMatrix <- matrix(0,nrow=numcat,ncol=numHrho)
    for(cc in 1:numHrho){
      transMatrix[which(unique1==cc),cc] <- 1
    }
    #minimal value for rho under H (without order constraints)
    pcat <- rep(1:length(ngroups),times=ngroups)
    LB <- unlist(unlist(lapply(1:numHrho,function(c){
      welk <- which(unique1==c)
      -1/(max(pvec[unlist(lapply(1:length(welk),function(j){
        which(pcat==welk[j])
      }))])-1)
    })))

    #initial value rho
    rhoH0ML <- unlist(lapply(1:numHrho,function(c){
      mean(rhoML[which(unique1==c)])
    }))
    rhoH <- rhoH0ML
    rhoS <- c(transMatrix%*%rhoH)
    RWsd <- rep(.1,numHrho)
    #initial unstandardized posterior (with uniform prior on rhoH)
    unpost <- int_lhood(rhoS,ngroups,pvec,N,K,Wmat1,zvec1,tWW2,tWz2,tzz2)

    #start sampling rho
    check1 <- 300 #check every 'check1' draws whether random walk sd needs to be increased/decreased
    acceptMat <- drawsMat <- matrix(0,nrow=samsize1,ncol=numHrho)
    #wer <- Sys.time()
    for(s in 1:samsize1){
      #sampling unique rho's under H
      for(j in 1:numHrho){
        #random walks
        rhoH_can <- rhoH
        rhoH_can[j] <- rnorm(1,mean=rhoH[j],sd=RWsd[j])
        if(rhoH_can[j]>LB[j] && rhoH_can[j]<1){
          #MH acceptance
          welk <- which(unique1==j)
          rhoS_can <- c(transMatrix%*%rhoH_can)
          unpost_can <- int_lhood(rhoS_can,ngroups,pvec,N,K,Wmat1,zvec1,tWW2,tWz2,tzz2)
          MHprob <- exp(unpost_can - unpost)
          if(runif(1) < MHprob){
            rhoH <- rhoH_can
            rhoS <- rhoS_can
            unpost <- unpost_can
            acceptMat[s,j] <- 1
          }
        }
      }
      drawsMat[s,] <- rhoH

      if(ceiling(s/check1)==s/check1){ #adjust sd of random walk based on acceptance proportions of last 'check1' draws
        probs <- apply(as.matrix(acceptMat[(s-check1+1):s,]),2,mean)
        upper1 <- .5
        lower1 <- .15
        RWsd[probs>upper1] <- RWsd[probs>upper1] * ( (probs[probs>upper1]-upper1)/(1-upper1) + 1)
        RWsd[probs<lower1] <- RWsd[probs<lower1] * 1/( 2 - (probs[probs<lower1])/lower1 )
      }
    }
    #Sys.time() - wer
    discard <- .2
    drawsMat <- as.matrix(drawsMat[(discard*samsize1+1):samsize1,]) #discard with 20% for burn-in
    samsize1 <- nrow(drawsMat)
    # get a tailored importance sampling estimate
    meanICC <- apply(drawsMat,2,mean)
    varICC <- apply(drawsMat,2,var)
    shape1IM <- - ((LB-meanICC)*((meanICC-1)*meanICC-meanICC*LB+LB+varICC)/((LB-1)*varICC))
    shape2IM <- - shape1IM*(meanICC-1)/(meanICC-LB)

    postestimates <- rbind(t(apply(drawsMat,2,mean)),
                           t(apply(drawsMat,2,median)),
                           apply(drawsMat,2,quantile,probs=c(.025,.975)))
    row.names(postestimates)[c(1,2)] <- c("mean","median")
    iccnames <- unlist(lapply(1:ncol(drawsMat),function(nc){paste0("icc",as.character(nc))}))
    colnames(postestimates) <- iccnames
    #compute posterior probs of positive icc of the free parameters
    postprobpositive <- apply(drawsMat>0,2,mean)
    priorprobpositive <- 1/(1 - LB)

    if(is.matrix(inequalities)){
      #compute prior probability of order constraints
      priordraws <- matrix(unlist(lapply(1:numHrho,function(cc){
        rbeta(samsize1,shape1=1,shape2=1)*(1-LB[cc])+LB[cc]
      })),ncol=numHrho)
      priorprob <- mean(apply(cbind(priordraws,rep(-1,samsize1))%*%t(inequalities) > 0,1,prod) )
      priorprob <- priorprob * (1-complement) + (1 - priorprob) * complement
      remove(priordraws)
      #compute posterior probability of order constraints
      postprob <- mean(apply(cbind(drawsMat,rep(-1,samsize1))%*%t(inequalities) > 0,1,prod) )
      postprob <- postprob * (1-complement) + (1 - postprob) * complement
    }else{
      priorprob <- postprob <- 1
    }
    #wer <- Sys.time()
    factor1 <- .6
    shape1IM <- shape1IM * factor1
    shape2IM <- shape2IM * factor1

    ISdraws <- matrix(unlist(lapply(1:numHrho,function(c){
      (rbeta(samsize2,shape1IM[c],shape2IM[c]) * (1 - LB[c]) + LB[c] ) * .99999
    })),ncol=numHrho)

    logintegrands <- unlist(lapply(1:samsize2,function(s){
      int_lhood(rhoS=c(transMatrix%*%ISdraws[s,]),ngroups,pvec,N,K,Wmat1,zvec1,tWW2,tWz2,tzz2) +
        sum(log(1/(1-LB))) -
        sum(dbeta((ISdraws[s,]-LB)/(1-LB),shape1=shape1IM,shape2=shape2IM,log=TRUE) +
              log(1/(1-LB)) )
    }))
    #Sys.time() - wer
    marglike <- log(mean(exp(logintegrands-max(logintegrands)))) +
      max(logintegrands) + log(postprob) - log(priorprob) + lgamma((N-K)/2) - (N-K)/2*log(pi)
  }
  return(list(marglike=marglike,postprob=postprob,priorprob=priorprob,
              postestimates=postestimates,postprobpositive=postprobpositive,
              priorprobpositive=priorprobpositive,postdraws=drawsMat))
}

Helmert = function(p){
  Helm <- diag(p)
  Helm[1,] <- 1/sqrt(p)
  for(pp in 2:p){
    Helm[pp,1:(pp-1)] <- 1/sqrt(pp*(pp-1))
    Helm[pp,pp] <- -(pp-1)/sqrt(pp*(pp-1))
  }
  return(Helm)
}

