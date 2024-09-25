### Joris Mulder 2019. Bayes factor for a one sample Student t test
### via adjusted FBFs (Mulder, 2014) using (m)lm-objects using a t.test object.

#' @method BF htest
#' @export
BF.htest <-
  function(x,
           hypothesis = NULL,
           prior.hyp.explo = NULL,
           prior.hyp.conf = NULL,
           prior.hyp = NULL,
           ...) {
    stop("Please use the function t_test() from the 'bain' package for a t-test or bartlett_test() from the 'BFpack' package for a test on group variances.")
  }



#' @importFrom stats approxfun
#' @describeIn BF BF S3 method for an object of class 't_test'
#' @method BF t_test
#' @export
BF.t_test <- function(x,
                      hypothesis = NULL,
                      prior.hyp.explo = NULL,
                      prior.hyp.conf = NULL,
                      prior.hyp = NULL,
                      complement = TRUE,
                      log = FALSE,
                      BF.type = NULL,
                      iter = 1e6,
                      ...){

  numpop <- length(x$estimate)

  if(is.null(BF.type)){
    BF.type <- "FBF"
  }
  if(!is.null(BF.type)){
    if(is.na(BF.type) | (BF.type!="FBF" & BF.type!="AFBF"))
      stop("The argument 'BF.type' must be 'FBF' (for the fractional BF) or 'AFBF' (for the adjusted fractional BF).")
  }
  if(BF.type=='AFBF'){
    bayesfactor <- "generalized adjusted fractional Bayes factors"
  }else{
    bayesfactor <- "generalized fractional Bayes factors"
  }

  logIN <- log

  if(numpop==1){ #one sample t test

    tvalue <- x$statistic
    mu0 <- x$null.value
    xbar <- x$estimate
    df <- x$parameter
    n <- df + 1
    stderr <- (xbar - mu0) / tvalue #standard error
    sigmaML <- stderr*sqrt(n-1)

    sampledata <- rnorm(n)
    sampledata <- (sampledata - mean(sampledata)) / sd(sampledata)
    sampledata <- sampledata * sqrt(x$v[1]) + xbar
    sampledata_explo <- sampledata - mu0
    mu <- rep(1,length(sampledata))
    df.draws <- data.frame(obs=sampledata,obs_explo=sampledata_explo,mu=mu)
    lm_conf <- lm(obs~mu-1,data=df.draws)
    BF_conf <- BF(lm_conf,
                  hypothesis=hypothesis,
                  prior.hyp.explo = prior.hyp.explo,
                  prior.hyp.conf = prior.hyp.conf,
                  prior.hyp=prior.hyp,
                  complement=complement,
                  BF.type=BF.type,
                  log=logIN)
    sampledata_explo <- sampledata - x$null.value
    lm_explo <- lm(obs_explo~mu-1,data=df.draws)
    BF_explo <- BF(lm_explo,
                   hypothesis=NULL,
                   prior.hyp.explo = prior.hyp.explo,
                   prior.hyp.conf = prior.hyp.conf,
                   prior.hyp=prior.hyp,
                   complement=complement,
                   BF.type=BF.type,
                   log=logIN)
    BF_out <- BF_conf
    BF_out$BFtu_exploratory <- BF_explo$BFtu_exploratory
    BF_out$PHP_exploratory <- BF_explo$PHP_exploratory
    colnames(BF_out$PHP_exploratory) <-
      c(paste0("Pr(=",as.character(mu0),")"),paste0("Pr(<",as.character(mu0),")"),
        paste0("Pr(>",as.character(mu0),")"))
    colnames(BF_out$BFtu_exploratory) <- c("BF0u","BF1u","BF2u")

  }else{ # two samples t test

    if(!grepl("Welch",x$method)){ # equal variances assumed
      # compute via a lm-analysis
      x1 <- c(rep(1,x$n[1]),rep(0,x$n[2]))
      y1 <- c(rep(0,x$n[1]),rep(1,x$n[2]))
      matx1y1 <- cbind(x1,y1)
      draw1 <- rnorm(x$n[1])
      out1 <- (draw1 - mean(draw1))/sd(draw1)*sqrt(x$v[1])+x$estimate[1]
      draw2 <- rnorm(x$n[2])
      out2 <- (draw2 - mean(draw2))/sd(draw2)*sqrt(x$v[2])+x$estimate[2]
      out <- c(out2,out1)
      # perform the test via a lm object using a factor for the group indicator
      # such that the name of the key variable (the difference between the means)
      # is called 'difference'
      df1 <- data.frame(out=out,differenc=factor(c(rep("a",x$n[2]),rep("e",x$n[1]))))
      lm1 <- lm(out ~ differenc,df1)
      BF_out <- BF(lm1,
                   hypothesis=hypothesis,
                   prior.hyp.explo=prior.hyp.explo,
                   prior.hyp.conf=prior.hyp.conf,
                   prior.hyp=prior.hyp,
                   complement=complement,
                   BF.type=BF.type,
                   log=logIN)
      BF_out$BFtu_exploratory <- t(as.matrix(BF_out$BFtu_exploratory[2,]))
      BF_out$PHP_exploratory <- t(as.matrix(BF_out$PHP_exploratory[2,]))
      row.names(BF_out$BFtu_exploratory) <- row.names(BF_out$PHP_exploratory) <- "difference"

      # if(!is.null(hypothesis)){
      #   BFtu_confirmatory <- BFlm1$BFtu_confirmatory
      #   PHP_confirmatory <- BFlm1$PHP_confirmatory
      #   BFmatrix_confirmatory <- BFlm1$BFmatrix_confirmatory
      #   BFtable <- BFlm1$BFtable_confirmatory
      #   hypotheses <- row.names(BFtable)
      #   priorprobs <- BFlm1$prior
      # }
    }else{ #equal variances not assumed. BF.lm cannot be used

      # check proper usage of argument 'prior.hyp.conf' and 'prior.hyp.explo'
      if(!is.null(prior.hyp.conf)){
        prior.hyp <- prior.hyp.conf
      }
      prior.hyp.explo <- process.prior.hyp.explo(prior_hyp_explo = prior.hyp.explo, model=x)

      diff.obs <- x$estimate[1] - x$estimate[2]
      nvec <- x$n
      bvec <- 2/nvec
      s2vec <- (x$v)*(x$n-1)
      #draws for Monte Carlo estimates
      sigma2_1.draws <- rgamma(iter,shape=(nvec[1]-1)/2,rate=s2vec[1]/2)
      sigma2b_1.draws <- rgamma(iter,shape=(nvec[1]*bvec[1]-1)/2,rate=s2vec[1]*bvec[1]/2)
      sigma2_2.draws <- rgamma(iter,shape=(nvec[2]-1)/2,rate=s2vec[2]/2)
      sigma2b_2.draws <- rgamma(iter,shape=(nvec[2]*bvec[2]-1)/2,rate=s2vec[2]*bvec[2]/2)
      relfit0 <- log(mean(dnorm(x$null.value,
                                mean=diff.obs,
                                sd=sqrt(1/sigma2_1.draws*1/nvec[1] + 1/sigma2_2.draws*1/nvec[2]))))
      relfit1 <- log(mean(pnorm(x$null.value,
                                mean=diff.obs,
                                sd=sqrt(1/sigma2_1.draws*1/nvec[1] + 1/sigma2_2.draws*1/nvec[2]))))
      relfit2 <- log(1 - exp(relfit1))
      if(BF.type == 'AFBF'){
        prior.mean <- x$null.value
      }else{
        prior.mean <- diff.obs
      }
      relcomp0 <- log(mean(dnorm(x$null.value,
                                 mean=prior.mean,
                                 sd=sqrt(1/sigma2b_1.draws*1/(nvec[1]*bvec[1]) + 1/sigma2b_2.draws*1/(nvec[2]*bvec[2])))))
      relcomp1 <- log(mean(pnorm(x$null.value,
                                 mean=prior.mean,
                                 sd=sqrt(1/sigma2b_1.draws*1/(nvec[1]*bvec[1]) + 1/sigma2b_2.draws*1/(nvec[2]*bvec[2])))))
      relcomp2 <- log(1 - exp(relcomp1))

      #exploratory Bayes factor test
      hypotheses_exploratory <- c("difference=0","difference<0","difference>0")
      logBFtu_exploratory <- c(relfit0-relcomp0,relfit1-relcomp1,relfit2-relcomp2)
      names(logBFtu_exploratory) <- hypotheses_exploratory
      BFtu_exploratory <- logBFtu_exploratory
      maxBFtu <- max(BFtu_exploratory)
      BFtu_explo_norm <- exp(BFtu_exploratory - maxBFtu) * prior.hyp.explo[[1]]
      PHP_exploratory <- matrix(BFtu_explo_norm/sum(BFtu_explo_norm),nrow=1)
      colnames(PHP_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
      row.names(PHP_exploratory) <- "difference"

      if(!is.null(hypothesis)){
        name1 <- "difference"
        parse_hyp <- parse_hypothesis(name1,hypothesis)
        parse_hyp$hyp_mat <- do.call(rbind, parse_hyp$hyp_mat)
        RrList <- make_RrList2(parse_hyp)
        RrE <- RrList[[1]]
        RrO <- RrList[[2]]

        if(length(RrE)==1){
          RrStack <- rbind(RrE[[1]],RrO[[1]])
          RrStack <- interval_RrStack(RrStack)
        }else{
          RrStack_list <- lapply(1:length(RrE),function(h){
            interval_RrStack(rbind(RrE[[h]],RrO[[h]]))
          })
          RrStack <- do.call(rbind,RrStack_list)
        }
        if(nrow(RrStack)>1){
          RStack <- RrStack[,-2]
          rStack <- RrStack[,2]
        }else{
          RStack <- matrix(RrStack[,-2],nrow=1)
          rStack <- RrStack[,2]
        }

        if(BF.type=='AFBF'){
          # check if a common boundary exists for prior location under all constrained hypotheses
          # necessary for the adjusted FBF
          if(nrow(RrStack) > 1){
            rref_ei <- rref(RrStack)
            nonzero <- rref_ei[,2]!=0
            if(max(nonzero)>0){
              row1 <- max(which(nonzero))
              if(sum(abs(rref_ei[row1,1]))==0){
                stop("No common boundary point for prior location. Conflicting constraints.")
              }
            }
          }
          prior.mean <- c(ginv(RStack)%*%rStack)
        }else{
          prior.mean <- diff.obs
        }

        relfit <- t(matrix(unlist(lapply(1:length(RrE),function(h){
          if(!is.null(RrE[[h]]) & is.null(RrO[[h]])){ #only an equality constraint
            nullvalue <- RrE[[h]][1,2]/RrE[[h]][1,1]
            relfit_h <- c(log(mean(dnorm(nullvalue,mean=diff.obs,
                                      sd=sqrt(1/sigma2_1.draws*1/nvec[1] + 1/sigma2_2.draws*1/nvec[2])))),
                          log(1))
            # relfit_h <- c(log(densN(nullvalue)),0)
          }else if(is.null(RrE[[h]]) & !is.null(RrO[[h]])){
            which.LB <- which(RrO[[h]][,1]>0)
            if(length(which.LB) == 0){
              LB <- -Inf
            }else{
              LB <- max(RrO[[h]][which.LB,2] / RrO[[h]][which.LB,1])
            }
            which.UB <- which(RrO[[h]][,1]<0)
            if(length(which.UB) == 0){
              UB <- Inf
            }else{
              UB <- min(RrO[[h]][which.UB,2] / RrO[[h]][which.UB,1])
            }
            relfit_h <- c(log(1),
                          log(mean(pnorm(UB,mean=diff.obs,
                                         sd=sqrt(1/sigma2_1.draws*1/nvec[1] + 1/sigma2_2.draws*1/nvec[2])) -
                                     pnorm(LB,mean=diff.obs,
                                           sd=sqrt(1/sigma2_1.draws*1/nvec[1] + 1/sigma2_2.draws*1/nvec[2]))
                                   )
                              )
                          )
          }else stop("hypothesis should either contain one equality constraint or inequality constraints on 'difference'.")
          return(relfit_h)
        })),nrow=2))
        #relfit <- exp(relfit)
        relcomp <- t(matrix(unlist(lapply(1:length(RrE),function(h){
          if(!is.null(RrE[[h]]) & is.null(RrO[[h]])){ #only an equality constraint
            nullvalue <- RrE[[h]][1,2]/RrE[[h]][1,1]
            relcomp_h <- c(log(mean(dnorm(nullvalue,mean=prior.mean,
                                          sd=sqrt(1/sigma2b_1.draws*1/(nvec[1]*bvec[1]) + 1/sigma2b_2.draws*1/(nvec[1]*bvec[1]))))),
                           log(1))
          }else if(is.null(RrE[[h]]) & !is.null(RrO[[h]])){ #order constraint(s)
            which.LB <- which(RrO[[h]][,1]>0)
            if(length(which.LB) == 0){
              LB <- -Inf
            }else{
              LB <- max(RrO[[h]][which.LB,2] / RrO[[h]][which.LB,1])
            }
            which.UB <- which(RrO[[h]][,1]<0)
            if(length(which.UB) == 0){
              UB <- Inf
            }else{
              UB <- min(RrO[[h]][which.UB,2] / RrO[[h]][which.UB,1])
            }
            relcomp_h <- c(log(1),
                          log(mean(pnorm(UB,mean=prior.mean,
                                         sd=sqrt(1/sigma2b_1.draws*1/(nvec[1]*bvec[1]) + 1/sigma2b_2.draws*1/(nvec[1]*bvec[1]))) -
                                     pnorm(LB,mean=prior.mean,
                                           sd=sqrt(1/sigma2b_1.draws*1/(nvec[1]*bvec[1]) + 1/sigma2b_2.draws*1/(nvec[1]*bvec[1])))
                                   )
                              )
                          )
          }else stop("hypothesis should either contain one equality constraint or inequality constraints on 'difference'.")
          return(relcomp_h)
        })),nrow=2))
        #relcomp <- exp(relcomp)
        row.names(relfit) <- row.names(relcomp) <- parse_hyp$original_hypothesis
        colnames(relfit) <- c("f=","f>")
        colnames(relcomp) <- c("c=","c>")
        if(complement == TRUE){
          #add complement to analysis
          if(sum(relcomp[,1]==0)>0){ #then there are order hypotheses
            # check if one-sided hypotheses are (partly) overlapping
            check.draws <- c(rnorm(iter,sd=1),rnorm(iter,sd=100))
            which.O <- which(unlist(lapply(RrO,function(x) {!is.null(x)})))
            checks <- Reduce("+",lapply(which.O,function(h){
              apply(as.matrix(RrO[[h]][,1]) %*% t(check.draws) > RrO[[h]][,2] %*%
                      t(rep(1,length(check.draws))),2,prod)
            }))
            if(max(checks)>1){ #(partly) overlapping
              warning("Complement hypothesis omitted because of (partly) overlapping hypotheses")
            }else{
              if(sum(exp(relcomp[which.O,2]))<1){#complement nonempty
                relcomp_c <- log(1-sum(exp(relcomp[which.O,2])))
                relcomp <- rbind(relcomp,c(0,relcomp_c))
                relfit_c <- log(1-sum(exp(relfit[which.O,2])))
                relfit <- rbind(relfit,c(0,relfit_c))
                row.names(relfit) <- row.names(relcomp) <- c(parse_hyp$original_hypothesis,"complement")
              }
            }
          }else{ #no order constraints
            relcomp <- rbind(relcomp,c(0,0))
            relfit <- rbind(relfit,c(0,0))
            row.names(relcomp) <- row.names(relfit) <- c(parse_hyp$original_hypothesis,"complement")
          }
        }
        # Check input of prior probabilies
        if(is.null(prior.hyp)){
          priorprobs <- rep(1/nrow(relcomp),nrow(relcomp))
        }else{
          if(!is.numeric(prior.hyp) || length(prior.hyp)!=nrow(relcomp)){
            warning(paste0("Argument 'prior.hyp' should be numeric and of length ",as.character(nrow(relcomp)),". Equal prior probabilities are used."))
            priorprobs <- rep(1/nrow(relcomp),nrow(relcomp))
          }else{
            priorprobs <- prior.hyp
          }
        }

        #compute Bayes factors and posterior probabilities for confirmatory test
        BFtu_confirmatory <- c(apply(relfit - relcomp, 1, sum))
        maxBFtu <- max(BFtu_confirmatory)
        PHP_confirmatory <- exp(BFtu_confirmatory-maxBFtu)*priorprobs / sum(exp(BFtu_confirmatory-maxBFtu)*priorprobs)
        BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)) -
          t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
        diag(BFmatrix_confirmatory) <- log(1)
        row.names(BFmatrix_confirmatory) <- colnames(BFmatrix_confirmatory) <- names(BFtu_confirmatory)
        relative_fit <- relfit
        relative_complexity <- relcomp

        BFtable <- cbind(relative_complexity,relative_fit,relative_fit[,1]-relative_complexity[,1],
                         relative_fit[,2]-relative_complexity[,2],apply(relative_fit,1,sum)/
                           apply(relative_complexity,1,prod),PHP_confirmatory)
        BFtable[,1:7] <- exp(BFtable[,1:7])
        row.names(BFtable) <- names(BFtu_confirmatory)
        colnames(BFtable) <- c("complex=","complex>","fit=","fit>","BF=","BF>","BF","PHP")
        hypotheses <- row.names(relative_complexity)

      }

      if(logIN == FALSE){
        BFtu_exploratory <- exp(BFtu_exploratory)
        if(!is.null(hypothesis)){
          BFtu_confirmatory <- exp(BFtu_confirmatory)
          BFmatrix_confirmatory <- exp(BFmatrix_confirmatory)
        }
      }

      if(is.null(hypothesis)){
        BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relative_fit <-
          relative_complexity <- BFtable <- hypotheses <- priorprobs <- NULL
      }

      rm(sigma2_1.draws);rm(sigma2_2.draws);rm(sigma2b_1.draws);rm(sigma2b_2.draws)
      parameter <- "mean difference"

      BF_out <- list(
        BFtu_exploratory=BFtu_exploratory,
        PHP_exploratory=PHP_exploratory,
        BFtu_confirmatory=BFtu_confirmatory,
        PHP_confirmatory=PHP_confirmatory,
        BFmatrix_confirmatory=BFmatrix_confirmatory,
        BFtable_confirmatory=BFtable,
        prior.hyp.explo=prior.hyp.explo,
        prior.hyp.conf=priorprobs,
        hypotheses=hypotheses,
        estimates=x$coefficients,
        model=x,
        bayesfactor=bayesfactor,
        parameter=parameter,
        log = logIN,
        fraction_number_groupIDs=2,
        fraction_groupID_observation=rep(1:2,times=nvec),
        call=match.call())
      class(BF_out) <- "BF"
    }
  }

  return(BF_out)
}













