### Joris Mulder 2019. Bayes factor testing of multiple random intercept models
### via multiple lmer-objects based on Mulder & Fox (2013, 2019).


#' @importFrom MCMCpack rinvgamma
# #' @importFrom lme4 getME
#' @method BF lmerMod
#' @export
BF.lmerMod <- function(x,
                   hypothesis = NULL,
                   prior = NULL,
                   parameter = NULL,
                   ...){
  # check if constrained hypotheses are formulated for confirmatory testing
  if(is.null(hypothesis)){
    constraints <- "exploratory"
  } else {
    constraints <- hypothesis
  }

  #prior probabilities of hypotheses
  if(is.null(prior)){
    priorprob <- "default"
  } else {
    priorprob <- prior
  }

  # ############ REPLACE THIS PART THAT ASSUMES AN INPUT OF A LIST OF LMERMOD-OBJECTS
  #
  # designX <- yobs <- group <- sorts <- list()
  # numcat <- length(x)
  # for(ca in 1:numcat){
  #   designX1 <- getME(x[[ca]],"X")
  #   yobs[[ca]] <- getME(x[[ca]],"y")
  #   group[[ca]] <- getME(x[[ca]],"flist")$Subject
  #   sorts <- order(group[[ca]])
  #   # sort data per category so that data from same groups are clustered over the rows
  #   yobs[[ca]] <- yobs[[ca]][sorts]
  #   designX1 <- designX1[sorts,]
  #   designX[[ca]] <- cbind(matrix(0,nrow=nrow(designX1),ncol=numcat),designX1[,-1])
  #   designX[[ca]][,ca] <- 1
  #   group[[ca]] <- group[[ca]][sorts]
  # }
  # ngroups <- rep(0,numcat)
  # for(ca in 1:numcat){
  #   # check that all groups/clusters are of equal size
  #   if(length(table(table(group[[ca]]))) != 1){stop("message")}
  #   # get number of groups per category
  #   ngroups[ca] <- length(levels(group[[ca]]))
  # }
  # for(ca in 1:(numcat-1)){
  #   if(names(table(table(group[[ca]]))) != names(table(table(group[[ca+1]])))){stop("message")}
  # }
  # p <- table(group[[ca]])[1]
  # names(p) <- NULL
  # ystack <- Reduce(c,yobs)
  # names(ystack) <- NULL
  # Xstack <- Reduce(rbind,designX)
  # row.names(Xstack) <- colnames(Xstack) <- NULL
  # yXstack <- cbind(ystack,Xstack)
  # #give names for icc's
  # iccnames <- unlist(lapply(1:numcat,function(nc){paste0("icc",as.character(nc))}))
  #
  # ############## END. REPLACE THIS PART

  # REPLACE THIS WITH PART WHERE DATA ARE NOT SORTED BY CLUSTERS
  # groups <- getME(x,"flist")[[1]]
  # if(length(table(table(groups)))>1){stop("Clusters are of unequal size.")}
  # p <- table(groups)[1]
  # # The data needs to be stacked per cluster.
  # ystack <- getME(x,"y")
  # Xstack <- getME(x,"X")
  # if(is.null(grouping_variables)){
  #   numcat <- 1
  #   ngroups <- nrow(Xstack)/p
  #   iccnames <- "icc"
  # }else{
  #   numcat <- length(grouping_variables)
  #   ngroups <- unlist(lapply(1:numcat,function(nc){
  #     sum(yX[,grouping_variables[nc]])/p
  #   }))
  #   iccnames <- unlist(lapply(1:numcat,function(nc){paste0("icc_",grouping_variables[nc])}))
  # }
  # END REPLACE 2


  #get names of the categories of clusters
  numcat <- length(x@cnms)
  namescat <- unlist(lapply(1:numcat,function(ca){
    x@cnms[[ca]]
  }))
  if(numcat>1){
    iccnames <- unlist(lapply(1:numcat,function(nc){paste0("icc_",namescat[nc])}))
  }else{ iccnames <- "icc" }

  # check if the lmer-model only has a random intercept or category specific random intercepts
  fixedeffectsonly <- F
  Zstack <- Reduce(cbind,getME(x,"mmList"))
  if(numcat>1){ #check if the random effects are category specific random intercepts
    for(ca in 1:numcat){
      freq_ca <- table(Zstack[,ca])
      if(sum(abs(sort(as.integer(names(freq_ca))) - c(0,1))) !=0){
        warning("only models with a single random intercept or category specific random intercepts are currently supported when testing icc's.")
        fixedeffectsonly <- T
      }
    }
    for(r in 1:nrow(Zstack)){
      if(table(Zstack[r,])["1"]!=1 || table(Zstack[r,])["0"]!=numcat-1){
        warning("only models with a single random intercept or category specific random intercepts are currently supported when testing icc's.")
        fixedeffectsonly <- T
      }
    }
  }else{
    freq_ca <- table(getME(x,"mmList")[[1]])
    if(as.integer(names(freq_ca))!=1 || length(as.integer(names(freq_ca)))!=1){
      warning("only models with a single random intercept or category specific random intercepts are currently supported when testing icc's.")
      fixedeffectsonly <- T
    }
  }

  # sort data per cluster
  clusterindex <- x@flist[[1]]
  if(length(table(table(clusterindex)))>1){stop("Clusters are of unequal size.")}
  p <- table(clusterindex)[1]
  levels(clusterindex) <- 1:length(levels(clusterindex))
  reorder1 <- order(as.integer(clusterindex))
  nclusters <- length(levels(clusterindex)) #total number of groups/clusters
  ystack <- getME(x,"y")[reorder1]
  Xstack <- getME(x,"X")[reorder1,]
  Zstack <- Zstack[reorder1,]
  # next sort data per category
  if(numcat>1){
    # sort rows per category
    catassign <- unlist(lapply((0:(nclusters-1))*p+1,function(cluster){
      which(Zstack[cluster,1:numcat]==1)
    }))
    ngroups <- table(catassign)
    names(ngroups) <- NULL
    reorder2 <- rep(unlist(lapply(1:numcat,function(ca){
      which(catassign==ca)
    }))-1,each=p)*p + rep(1:p,nclusters)
    Xstack <- Xstack[reorder2,]
    ystack <- ystack[reorder2]
    Zstack <- Zstack[reorder2,]
  }else{ ngroups <- nclusters }

  # marginal likelihoods of unconstrained model and unconstrained estimation
  #default prior for icc's is stretched beta(1,1)
  if(fixedeffectsonly==F){
    shape0 <- c(1,1)
    cat("First, unconstrained icc analysis...")
    cat("\n")
    cat("\n")
    marglike_Hu <- marglike2_Hq(cbind(ystack,Xstack),ngroups,p,shape1=shape0[1],shape2=shape0[2],
                                         samsize1=5e3,samsize2=5e3,unique=1:numcat,inequalities=0)
    postestimates <- marglike_Hu[[4]]
    colnames(postestimates) <- iccnames

    # exploratory testing
    cat("Next, Bayes factor computation for exploratory testing of icc's...")
    cat("\n")
    BFtu_exploratory_icc <- t(matrix(unlist(lapply(1:numcat,function(nc){

      cat(paste0(iccnames[nc],"; "))
      cat("\n")

      if(numcat>1){
        unique_c <- rep(1,numcat)
        unique_c[nc] <- 0
        unique_c[-nc] <- 1:(numcat-1)
      }else unique_c <- 0

      marglike_explo <- rep(0,3)
      # zero icc
      marglike_explo[1] <- marglike2_Hq(cbind(ystack,Xstack),ngroups,p,shape1=shape0[1],
                                                 shape2=shape0[2],unique=unique_c)[[1]]
      inequalities2 = matrix(c(unique_c==0,0),ncol=numcat+1)
      # positive icc
      marglike_positive <- marglike2_Hq(cbind(ystack,Xstack),ngroups,p,shape1=shape0[1],
                                                 shape2=shape0[2],unique=1:numcat,inequalities=inequalities2)
      marglike_explo[3] <- marglike_positive[[1]]
      # negative icc
      marglike_explo[2] <- marglike_positive[[1]] - log(marglike_positive[[2]]) + log(marglike_positive[[3]]) +
        log(1-marglike_positive[[2]]) - log(1-marglike_positive[[3]])
      return(exp(marglike_explo - marglike_Hu[[1]]))
    })),nrow=3))
    colnames(BFtu_exploratory_icc) <- c("icc=0","icc<0","icc>0")
    row.names(BFtu_exploratory_icc) <- iccnames
    PHP_exploratory_icc <- round(BFtu_exploratory_icc / apply(BFtu_exploratory_icc,1,sum),3)
    priorprobs <- rep(1,3)/3 #prior probs for exploratory tests
    cat("\n")
    cat("\n")
    if(constraints!="exploratory"){ # confirmatory test with constrained hypotheses on icc's.
      cat("Bayes factor computation for confirmatory testing of icc's...")
      cat("\n")
      parse_hyp <- parse_hypothesis(iccnames,constraints)
      RrList <- make_RrList2(parse_hyp)
      RrE <- RrList[[1]]
      RrO <- RrList[[2]]
      # check if icc's are only tested against each other or against zero
      numhyp <- length(RrE)
      for(h in 1:numhyp){
        if(!is.null(RrE[[h]])){
          for(r in 1:nrow(RrE[[h]])){
            row1 <- RrE[[h]][r,]
            if( !(sum(row1)==1 || sum(row1)==0) ){
              stop("icc's can only be compared with each other or to zero.")
            }
          }
        }
        if(!is.null(RrO[[h]])){
          for(r in 1:nrow(RrO[[h]])){
            freq1 <- table(sort(RrO[[h]][r,]))
            row1 <- RrO[[h]][r,]
            if( !(sum(row1)==1 || sum(row1)==0) ){
              stop("icc's can only be compared with each other or to zero.")
            }
          }
        }
      }

      output_marglike_icc <- t(matrix(unlist(lapply(1:numhyp, function(h){

        cat(paste0(parse_hyp$original_hypothesis[h],"; "))
        cat("\n")

        # code equal icc's with same integer for marglike2_Hq function
        unique_h <- 1:numcat
        if(!is.null(RrE[[h]])){
          for(r in 1:nrow(RrE[[h]])){
            row1 <- RrE[[h]][r,]
            if(sum(row1)==0){
              which_equal <- which(row1!=0)
              which_highercodegroup <- which(unique_h > max(unique_h[which_equal]))
              unique_h[which_equal] <- min(unique_h[which_equal])
              unique_h[which_highercodegroup] <- unique_h[which_highercodegroup] - 1
            } else if(sum(row1)==1){
              which_zero <- which(row1==1)
              which_highercodegroup <- which(unique_h > max(unique_h[which_zero]))
              unique_h[which_zero] <- 0
              unique_h[which_highercodegroup] <- unique_h[which_highercodegroup] - 1
            }
          }
        }
        if(!is.null(RrO[[h]])){
          inequalities_h <- matrix(0,nrow(RrO[[h]]),ncol=length(unique(unique_h))+1)
          for(u in sort(unique(unique_h[unique_h>0]))){
            inequalities_h[,u] <- apply(as.matrix(RrO[[h]][,which(unique_h == u)]),1,sum)
          }
        } else inequalities_h = 0

        marglike2_h <- marglike2_Hq(cbind(ystack,Xstack),ngroups,p,shape1=shape0[1],shape2=shape0[2],samsize1=5e3,samsize2=5e3,
                                    unique=unique_h,inequalities=inequalities_h)[1:3]
        return(c(unlist(marglike2_h),ifelse(is.null(RrE[[h]]),1,0)))
      })),nrow=4))
      if(sum(output_marglike_icc[,4])==0){ #the complement is equivalent to the unconstrained model
        output_marglike_icc <- rbind(output_marglike_icc,c(unlist(marglike_Hu),1))
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
      BFtu_confirmatory_icc <- exp(output_marglike_icc[,1] - marglike_Hu[[1]])

      #compute BFmatrix and PHPs
      logBFmatrix <- matrix(rep(output_marglike_icc[,1],numhyp+1),nrow=numhyp+1) -
        matrix(rep(output_marglike_icc[,1],each=numhyp+1),nrow=numhyp+1)
      row.names(logBFmatrix) <- colnames(logBFmatrix) <- c(parse_hyp$original_hypothesis,"complement")
      BFmatrix_confirmatory_icc <- round(exp(logBFmatrix),3)
      BFta_confirmatory_icc <- exp(output_marglike_icc[,1] - max(output_marglike_icc[,1]))
      # Change prior probs in case of default setting
      if(priorprob=="default"){
        priorprobs <- rep(1/length(BFta_confirmatory_icc),length(BFta_confirmatory_icc))
      }else{
        priorprobs <- priorprobs/sum(priorprobs)
      }
      PHP_confirmatory_icc <- round(priorprobs*BFta_confirmatory_icc / sum(priorprobs*BFta_confirmatory_icc),3)
    }else{
      BFmatrix_confirmatory_icc <- PHP_confirmatory_icc <- BFtu_confirmatory_icc <- NULL
    }
  }else{
    BFtu_exploratory_icc <- PHP_exploratory_icc <- BFtu_confirmatory_icc <-
      PHP_confirmatory_icc <- BFmatrix_confirmatory_icc <- postestimates <-
      priorprobs <- NULL
  }

  #### Code for testing the fixed effects here.

  #### End code for testing the fixed effects here.

  #####
  # Combine results of tests of fixed effects and
  #####

  BFlm_out <- list(
    BFtu_exploratory=BFtu_exploratory_icc,
    PHP_exploratory=PHP_exploratory_icc,
    BFtu_confirmatory=BFtu_confirmatory_icc,
    PHP_confirmatory=PHP_confirmatory_icc,
    BFmatrix_confirmatory=BFmatrix_confirmatory_icc,
    relative_fit=NULL,
    relative_complexity=NULL,
    model=x,
    ngroups=ngroups,
    p=p,
    Xstack=Xstack,
    Zstack=Zstack,
    ystack=ystack,
    constraints=constraints,
    priorprobs=priorprobs,
    estimates=postestimates,
    categories=namescat)

  class(BFlm_out) <- "BF"

  return(BFlm_out)
}

BFupdatelmer <- function(x,
                         ...){
  stop("REMINDER: still needs to be finished.")
}


# Functions that are called when computing marginal likelihoods for constrained
# hypotheses on icc's / between-groups variances

#create Helmert matrix
Helmert = function(p){
  Helm <- diag(p)
  Helm[1,] <- 1/sqrt(p)
  for(pp in 2:p){
    Helm[pp,1:(pp-1)] <- 1/sqrt(pp*(pp-1))
    Helm[pp,pp] <- -(pp-1)/sqrt(pp*(pp-1))
  }
  return(Helm)
}
#sampler under a equality-constrained hypothesis on icc's
Gibbs2 <- function(zW,ngroups,p,shape1,shape2,bB,bW,unique,T0,V1,inequalities=0,samsize=3e3){
  #Gibbs sampler of model parameters under the constrained hypothesis H*_q,
  #excluding the inequality constraints.
  #returns hyper parameters of a shifted-beta(shape1,shape2,-1/(p-1),1)-
  #distribution used for importance sampling.

  clusters <- length(ngroups)
  N <- sum(ngroups)
  K <- ncol(zW)-1

  Hp <- Helmert(p)
  select1 <- p*(0:(N-1))+1

  Wmat <- matrix(zW[,-1],ncol=ncol(zW)-1)
  zvec <- zW[,1]

  #initial values
  sigma2 <- 1
  tauV <- rep(.5,V1)
  beta <- rep(1,K)
  psi <- rep(1,V1)

  transMatrix <- matrix(0,ncol=V1,nrow=clusters)
  if(T0<clusters){
    for(cc in (T0+1):clusters){
      transMatrix[cc,unique[cc]] <- 1
    }
  }

  T0check <- T0>0
  psi.check <- shape1>0
  psi <- psi * psi.check

  #store ICC's
  rhoMat2 <- matrix(0,ncol=V1,nrow=samsize)

  burnin <- 1e3
  #  pb = txtProgressBar(min = 0, max = burnin, initial = 0)

  #  print("burn-in")

  #burnin
  for(ss in 1:burnin){
    #1. draw beta | sigma2, tau, psi, y
    tauC <- c(transMatrix%*%tauV)
    vars <- unlist(lapply(1:clusters,function(cc){
      rep(c((sigma2+p*tauC[cc])*bB[cc]**(-1),rep(sigma2*bW**(-1),p-1)),ngroups[cc])
    }))
    covBeta <- solve(t(Wmat/vars)%*%Wmat) #+ diag(ncol(Wmat))*.00001
    meanBeta <- covBeta%*%t(Wmat/vars)%*%zvec
    beta <- c(rmvnorm(1,mean=meanBeta,sigma=covBeta))
    #beta = c(0,0)

    #2. draw (sigma2,psi,tau) | beta, y
    diffs <- zvec - Wmat%*%beta
    sumsquares.tau <- unlist(lapply(1:clusters,function(cc){
      select_c <- p*((sum(ngroups[1:cc])-ngroups[cc]):(sum(ngroups[1:cc])-1))+1
      sum(diffs[select_c]**2)
    }))
    #2a. draw sigma2| beta, y
    scale.sigma2 <- sum((diffs[-select1])**2)*bW/2 +
      sum(sumsquares.tau[1:T0]*bB[1:T0])/2*T0check
    shape.sigma2 <- bW*N*(p-1)/2 + sum(bB[1:T0]*ngroups[1:T0])/2*T0check
    sigma2 <- MCMCpack::rinvgamma(1,shape=shape.sigma2,scale=scale.sigma2)

    #2b. draw psi | sigma2, beta, y
    # if psi.check==F then psi is set to zero
    shape.psi <- shape1*psi.check + (1-psi.check)*rep(10,V1)
    rate.psi <- p/((p-1)*sigma2)
    psi <- rgamma(V1,shape=shape.psi,rate=rate.psi)*psi.check

    #2c. draw tau | sigma2, psi, beta, y
    scale.tau <- c(t(sumsquares.tau*bB)%*%transMatrix)/(2*p) + psi
    shape.tau = c(t(bB*ngroups)%*%transMatrix)/2 + shape2
    tauV <- MCMCpack::rinvgamma(V1,shape=shape.tau,scale=scale.tau) - sigma2/p #+ .00001

    #    setTxtProgressBar(pb,ss)
  }

  #  cat("\n")
  #  print("finding good proposal density")

  #  pb = txtProgressBar(min = 0, max = samsize, initial = 0)

  for(ss in 1:samsize){
    #1. draw beta | sigma2, tau, psi, y
    tauC <- c(transMatrix%*%tauV)
    vars <- unlist(lapply(1:clusters,function(cc){
      rep(c((sigma2+p*tauC[cc])*bB[cc]**(-1),rep(sigma2*bW**(-1),p-1)),ngroups[cc])
    }))
    covBeta <- solve(t(Wmat/vars)%*%Wmat) #+ diag(ncol(Wmat))*.00001
    meanBeta <- covBeta%*%t(Wmat/vars)%*%zvec
    beta <- c(rmvnorm(1,mean=meanBeta,sigma=covBeta))

    #2. draw (sigma2,psi,tau) | beta, y
    diffs <- c(zvec - Wmat%*%beta)
    sumsquares.tau <- unlist(lapply(1:clusters,function(cc){
      select_c <- p*((sum(ngroups[1:cc])-ngroups[cc]):(sum(ngroups[1:cc])-1))+1
      sum(diffs[select_c]**2)
    }))
    #2a. draw sigma2| beta, y
    scale.sigma2 <- sum((diffs[-select1])**2)*bW/2 +
      sum(sumsquares.tau[1:T0]*bB[1:T0])/2*T0check
    shape.sigma2 <- bW*N*(p-1)/2 + sum(bB[1:T0]*ngroups[1:T0])/2*T0check
    sigma2 <- MCMCpack::rinvgamma(1,shape=shape.sigma2,scale=scale.sigma2)

    #2b. draw psi | sigma2, beta, y
    # if psi.check==F then psi is set to zero
    shape.psi <- shape1*psi.check + (1-psi.check)*rep(10,V1)
    rate.psi <- p/((p-1)*sigma2)
    psi <- rgamma(V1,shape=shape.psi,rate=rate.psi)*psi.check

    #2c. draw tau | sigma2, psi, beta, y
    scale.tau <- c(t(sumsquares.tau*bB)%*%transMatrix)/(2*p) + psi
    shape.tau <- c(t(bB*ngroups)%*%transMatrix)/2 + shape2
    tauV <- MCMCpack::rinvgamma(V1,shape=shape.tau,scale=scale.tau) - sigma2/p #+ .00001

    #2d. compute rho
    rho <- tauV/(tauV+sigma2)

    rhoMat2[ss,] <- rho
  }

  meanICC <- apply(rhoMat2,2,mean)
  varICC <- apply(rhoMat2,2,var)
  LB <- -1/(p-1)
  shape1IM <- (meanICC*(p-1)+1)*(meanICC**2*(1-p)+meanICC*(p-2)-varICC*(p-1)+1)/((p-1)*p*varICC)
  shape2IM <- (meanICC-1)*(meanICC^2*(p-1)-meanICC*(p-2)+varICC*(p-1)-1)/(varICC*p)

  post.prob <- 1
  if(sum(abs(inequalities))>0){
    inequalities0 <- matrix(inequalities[,1:V1],ncol=V1)
    post.prob <- mean(apply(rhoMat2%*%t(inequalities0)>rep(1,samsize)%*%t(inequalities[,V1+1]),1,prod))
  }

  return(list(shape1IM,shape2IM,post.prob,rhoMat2))
}
#integrand of the marginal likelhood
logintegrand_Hq <- function(rhoV,zvec,Wmat,p,ngroups,shape1,shape2,bB,bW,transMatrix){
  #zvec is the stacked Helmert transformed outcome vector.
  #Wmat is the stacked Helmert transformed matrix with predictors.
  ngroupsVb <- c((ngroups*bB)%*%transMatrix)
  rhoC <- c(transMatrix%*%rhoV)
  N <- sum(ngroups)

  clusters <- length(ngroups)
  K <- ncol(Wmat)

  #the normalizing constants (term1 and term2) will be included after averaging.
  term3 <- sum( (-ngroupsVb/2+shape1-1)*log(1+(p-1)*rhoV) +
                 (ngroupsVb/2+shape2-1)*log(1-rhoV) )

  vars <- unlist(lapply(1:clusters,function(cc){
    rep( c( (1+(p-1)*rhoC[cc])/(1-rhoC[cc])*bB[cc]**(-1),rep(bW**(-1),p-1)),ngroups[cc]) + .0000001
  }))
  tWDi <- t(Wmat/vars)
  tWDiW <- tWDi%*%Wmat
  tWDiWi <- solve(tWDiW)
  betaHat <- c(tWDiWi%*%tWDi%*%zvec)
  diffHat <- c(zvec - Wmat%*%betaHat)
  s2Hat <- sum(diffHat*diffHat/vars)

  term4 <- -.5*(sum(ngroups*bB)+N*(p-1)*bW-K)*log(s2Hat)
  term5 <- -.5*log(det(tWDiW))

  return(#term1+term2+
    term3 + term4 + term5)
}
#log density of stretched-beta distribution in (-1/(p-1);1)
log_dshiftedbeta <- function(x,shape1,shape2,p){
  lgamma(shape1+shape2)-lgamma(shape1)-lgamma(shape2)+(-shape1-shape2+1)*log(p)+shape2*log(p-1)+(shape1-1)*log(1+(p-1)*x)+(shape2-1)*log(1-x)
}
#computation of the marginal likelihood of Hq.
marglike_Hq <- function(yX,ngroups,p,shape1=1,shape2=1,samsize1=5e3,samsize2=5e3,bB=1,bW=1,unique,inequalities=0,complement=FALSE){

  #E.g., for clusters=5, unique=c(0,2,0,1,1),inequalities=[1 -1 0],
  #the hypothesis equals, Hq:rho1=rho3=0, rho4=rho5>rho2
  #samsize1 sets the sample size for the number of draws from the proposal distribution for the importance
  #sample estimate of the marginal likelihood excluding the inequality constraints.
  #samsize2 sets the sample size for computing the probability that the inequality constraints hold,
  #and for constructing the proposal distribution.

  clusters <- length(ngroups)
  N <- sum(ngroups)
  Hp <- Helmert(p)

  if(length(bB)==1){
    bB <- rep(bB,clusters)
  }

  #re-order: first zero-rho's, then rho's coded as 1, then 2, etc.
  T0 <- sum(unique==0) #number of zero tau's
  V1 <- max(unique)
  K1 <- ncol(yX)-1
  lvec <- rep(0,V1)
  kvec <- rep(0,V1)
  ord1 <- rep(0,clusters)
  if(T0>0){
    ord1[1:T0] <- which(unique==0)
  }

  plek <- T0
  for(hh in 1:V1){
    welke_h <- which(unique==hh)
    lvec[hh] <- length(welke_h)
    kvec[hh] <- sum(unique<hh)
    ord1[plek+1:length(welke_h)] <- welke_h
    plek <- plek + length(welke_h)
  }
  #Change order of the data
  yXdummy <- yX
  location <- 0
  for(cc in 1:clusters){
    welke <- ord1[cc]
    maxcc <- sum(ngroups[1:welke])*p
    mincc <- (sum(ngroups[1:welke])-ngroups[welke])*p+1
    yXdummy[location+1:(maxcc-mincc+1),] <- yX[mincc:maxcc,]
    location <- location + maxcc-mincc+1
  }
  yX <- yXdummy
  ngroups <- ngroups[ord1]
  bB <- bB[ord1]
  unique <- unique[ord1]
  T0check <- (T0>0)

  if(length(shape1)==1){
    shape1 <- rep(shape1,V1)
  }
  if(length(shape2)==1){
    shape2 <- rep(shape2,V1)
  }

  transMatrix <- matrix(0,ncol=V1,nrow=clusters)
  if(T0<clusters){
    for(cc in (T0+1):clusters){
      transMatrix[cc,unique[cc]] <- 1
    }
  }

  zW <- yX
  for(ii in 1:N){
    zW[((ii-1)*p+1):(ii*p),] <- Hp%*%yX[((ii-1)*p+1):(ii*p),]
  }
  Xmat <- matrix(yX[,-1],ncol=ncol(yX)-1)
  yvec <- yX[,1]
  Wmat <- matrix(zW[,-1],ncol=ncol(zW)-1)
  zvec <- zW[,1]

  inequalities0 <- 0
  if(sum(abs(inequalities))>0){
    inequalities0 <- matrix(inequalities[,1:V1],ncol=V1)
  }

  #determine hyperparameters for the importance sampler and prob of inequality constraints Hq
  out2 <- Gibbs2(zW,ngroups,p,shape1,shape2,bB,bW,unique,T0,V1,inequalities,samsize=samsize2)
  post.prob.Hq <- out2[[3]]*(1-complement) + (1-out2[[3]])*complement
  # estimates under the marginal model retricted with the equality constraints
  postestimates <- rbind(t(apply(out2[[4]],2,mean)),
                         t(apply(out2[[4]],2,median)),
                         apply(out2[[4]],2,quantile,probs=c(.025,.975)))
  row.names(postestimates)[c(1,2)] <- c("mean","median")
  iccnames <- unlist(lapply(1:ncol(out2[[4]]),function(nc){paste0("icc",as.character(nc))}))
  colnames(postestimates) <- iccnames

  #compute the normalizing constant in the prior
  logKq <- 0
  prior.prob.Hq <- 1
  LB <- -1/(p-1)
  if(prod(shape1*shape2)>0){
    logKq <- sum(lgamma(shape1+shape2)-lgamma(shape1)-lgamma(shape2)) +
      sum(shape2*log(p-1) - (shape1+shape2-1)*log(p))
    if(sum(abs(inequalities))>0){
      priordraws <- matrix(unlist(lapply(1:V1,function(cc){
        rbeta(samsize2,shape1=shape1[cc],shape2=shape2[cc])*(1-LB)+LB
      })),ncol=V1)
      prior.prob.Hq <- mean(apply(priordraws%*%t(inequalities0)>rep(1,samsize2)%*%t(inequalities[,V1+1]),1,prod))
      prior.prob.Hq <- prior.prob.Hq*(1-complement) + (1-prior.prob.Hq)*complement
    }
  }

  #multiply with .6 to make the importance sampler distribution slightly wider
  #draws from proposal density
  factor1 <- .6
  hyperIS <- list(out2[[1]]*factor1,out2[[2]]*factor1)
  rdraws <- samsize1
  LB <- -1/(p-1)
  rhodraws <- matrix(unlist(lapply(1:V1,function(cc){
    (rbeta(rdraws,shape1=hyperIS[[1]][cc],shape2=hyperIS[[2]][cc])*(1-LB)+LB)*.99999
  })),ncol=V1)
  #multiplication with .9999 to avoid points too close to boundary

  #compute prior*likelihood for importance sample draws
  logintegrands <- unlist(lapply(1:rdraws,function(ss){
    logintegrand_Hq(rhoV=rhodraws[ss,],zvec,Wmat,p,ngroups,shape1,shape2,bB,bW,transMatrix) -
      sum(log_dshiftedbeta(rhodraws[ss,],shape1=hyperIS[[1]],shape2=hyperIS[[2]],p))
  }))

  term1 <- -.5*(sum(ngroups*bB)+N*(p-1)*bW-K1)*log(pi)
  term2 <- lgamma(.5*(sum(ngroups*bB)+N*(p-1)*bW-K1))
  logintegrands <- logintegrands + term1 + term2 + logKq - log(prior.prob.Hq)
  marglike.Hq <- log(mean(exp(logintegrands-max(logintegrands)))) +
    max(logintegrands) + log(post.prob.Hq)

  return(list(marglike.Hq,post.prob.Hq,prior.prob.Hq,postestimates))
}
#marginal likelihood for H0:rho1=...=rhoC=0
marglike_H0 <- function(yX,ngroups,p,bB=1,bW=1){

  N <- sum(ngroups)
  clusters <- length(ngroups)
  K <- ncol(yX)-1
  Hp <- Helmert(p)

  if(length(bB)==1){
    bB <- rep(bB,clusters)
  }

  zW <- yX
  for(ii in 1:N){
    zW[((ii-1)*p+1):(ii*p),] = Hp%*%yX[((ii-1)*p+1):(ii*p),]
  }
  Xmat <- yX[,-1]
  yvec <- yX[,1]
  Wmat <- zW[,-1]
  zvec <- zW[,1]

  term1 <- -.5*(sum(ngroups*bB)+N*(p-1)*bW-K)*log(pi)
  term2 <- lgamma(.5*(sum(ngroups*bB)+N*(p-1)*bW-K))

  vars <- unlist(lapply(1:clusters,function(cc){
    rep(c(bB[cc]**(-1),rep(bW**(-1),p-1)),ngroups[cc])
  }))

  tWDi <- t(Wmat/vars)
  tWDiW <- tWDi%*%Wmat
  tWDiWi <- solve(tWDiW)
  betaHat <- c(tWDiWi%*%tWDi%*%zvec)
  diffHat <- c(zvec - Wmat%*%betaHat)
  s2Hat <- sum(diffHat*diffHat/vars)

  term3 <- -.5*(sum(ngroups*bB)+N*(p-1)*bW-K)*log(s2Hat)
  term4 <- -.5*log(det(tWDiW))

  return(list(term1+term2+term3+term4,1,1,NA))
}
#wrapper for marglike_Hq but applies FBF approach with minimal fractions if shape1=shape2=0
marglike2_Hq <- function(yX,ngroups,p,shape1=1,shape2=1,samsize1=5e3,samsize2=5e3,unique,inequalities=0,complement=FALSE){
  if(sum(abs(unique))>0){
    if(prod(shape1*shape2)==0){#apply FBF methodology
      #default fractions need to be checked...
      bBmin <- 2/ngroups # to identify the fixed cluster specific intercepts and cluster specific tau's
      bWmin <- (ncol(yX)-length(ngroups))/(sum(ngroups)*(p-1)) # to identify sigma2 and the remaining fixed effects
      pyb <- marglike_Hq(yX,ngroups,p,shape1,shape2,samsize1,samsize2,bB=bBmin,bW=bWmin,unique,inequalities,complement)
      py1 <- marglike_Hq(yX,ngroups,p,shape1,shape2,samsize1,samsize2,bB=1,bW=1,unique,inequalities,complement)
      outpHq <- py1[[1]] - pyb[[1]]
    }else{
      outpHq <- marglike_Hq(yX,ngroups,p,shape1,shape2,samsize1,samsize2,bB=1,bW=1,unique,inequalities,complement)
    }
  }else{
    if(prod(shape1*shape2)==0){#apply FBF methodology
      #default fractions need to be checked...
      bBmin <- 2/ngroups # to identify the fixed cluster specific intercepts and cluster specific tau's
      bWmin <- (ncol(yX)-length(ngroups))/(sum(ngroups)*(p-1)) # to identify sigma2 and the remaining fixed effects
      pyb <- marglike_H0(yX,ngroups,p,bB=bBmin,bW=bWmin)
      py1 <- marglike_H0(yX,ngroups,p,bB=1,bW=1)
      outpHq <- py1[[1]] - pyb[[1]]
    }else{
      outpHq <- marglike_H0(yX,ngroups,p,bB=1,bW=1)
    }
  }
  return(outpHq)
}
