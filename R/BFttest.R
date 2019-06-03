### Joris Mulder 2019. Bayes factor for a one sample Student t test
### via adjusted FBFs (Mulder, 2014) using (m)lm-objects using a t.test object.

#' @method bain htest
#' @export
bain.htest <-
  function(x,
           hypothesis,
           prior = NULL,
           parameter = NULL,
           ...) {
    stop("The standard t.test() function from the 'stats' package does not return variance and sample size, which are required to run bain. Please use the function t_test() from the 'bain' package instead. It accepts the same arguments.")
  }

BF.bain_htest <- function(x,
                      hypothesis = NULL,
                      prior = NULL,
                      parameter = NULL,
                      ...){

  numpop <- length(x$estimate)

  if(numpop==1){ #one sample t test
    # tvalue <- x$statistic
    # names(tvalue) <- NULL
    # n <- x$parameter + 1
    # names(n) <- NULL
    # mu0 <- x$null.value
    # names(mu0) <- NULL
    # xbar <- x$estimate
    # names(xbar) <- NULL
    # sd1 <- sqrt(n) * (xbar - mu0) / tvalue
    # y1 <- rnorm(n)
    # y1 <- y1 - mean(y1)
    # y1 <- sd1*y1/sd(y1) + xbar - mu0
    # lm1 <- lm(y1 ~ 1)
    # names(lm1$coefficients) <- "mu_minus_mu0"
    #
    #
    tvalue <- x$statistic
    mu0 <- x$null.value
    xbar <- x$estimate
    df <- x$parameter
    n <- df + 1
    stderr <- (xbar - mu0) / tvalue #standard error
    sigmaML <- stderr*sqrt(n-1)
    relfit0 <- dt((xbar-mu0)/stderr,df=df,log=TRUE) - log(stderr)
    relfit1 <- log(1-pt((xbar-mu0)/stderr,df=df))
    relfit2 <- log(pt((xbar-mu0)/stderr,df=df))
    relcomp0 <- dt(0,df=1,log=T) - log(sigmaML)
    relcomp1 <- relcomp2 <- log(.5)

    if(is.null(hypothesis)){
      if(x$alternative=="two.sided"){
        if(!grepl("Paired",x$method)){
          hypotheses <- c(paste0("true mean is equal to ",as.character(mu0)),paste0("true mean is not equal to ",as.character(mu0)))
        }else{
          hypotheses <- c("true difference is equal to 0","true difference is not equal to 0")
        }
        logBFtu <- c(relfit0-relcomp0,0)
        names(logBFtu) <- hypotheses
        BFtu <- exp(logBFtu)
        PHP <- BFtu/sum(BFtu)
        BFmatrix <- matrix(rep(BFtu,2),ncol=2) / matrix(rep(BFtu,2),ncol=2,byrow=T)
        logBFmatrix <- matrix(rep(logBFtu,2),ncol=2) - matrix(rep(logBFtu,2),ncol=2,byrow=T)
        row.names(BFmatrix) <- row.names(logBFmatrix) <- hypotheses
        colnames(BFmatrix) <- colnames(logBFmatrix) <- hypotheses
        relfit <- matrix(c(exp(relfit0),rep(1,3)),ncol=2)
        relcomp <- matrix(c(exp(relcomp0),rep(1,3)),ncol=2)
        row.names(relfit) <- row.names(relcomp) <- hypotheses
        colnames(relfit) <- c("f=","f>")
        colnames(relcomp) <- c("c=","c>")

        # BF1 <- BF(lm1,hypothesis=paste0("mu_minus_mu0=",as.character(mu0)),
        #           prior=prior)
      } else if(x$alternative=="less"){
        if(!grepl("Paired",x$method)){
          hypotheses <- c(paste0("true mean is greater than or equal to ",as.character(mu0)),paste0("true mean is less than ",as.character(mu0)))
        }else{
          hypotheses <- c("true difference is greater than or equal to 0","true difference is less than 0")
        }
        logBFtu <- c(relfit2-relcomp2,relfit1-relcomp1)
        names(logBFtu) <- hypotheses
        BFtu <- exp(logBFtu)
        PHP <- BFtu/sum(BFtu)
        BFmatrix <- matrix(rep(BFtu,2),ncol=2) / matrix(rep(BFtu,2),ncol=2,byrow=T)
        logBFmatrix <- matrix(rep(logBFtu,2),ncol=2) - matrix(rep(logBFtu,2),ncol=2,byrow=T)
        row.names(BFmatrix) <- row.names(logBFmatrix) <- hypotheses
        colnames(BFmatrix) <- colnames(logBFmatrix) <- hypotheses
        relfit <- matrix(c(rep(1,2),exp(relfit2),exp(relfit1)),ncol=2)
        relcomp <- matrix(c(rep(1,2),rep(.5,2)),ncol=2)
        row.names(relfit) <- row.names(relcomp) <- hypotheses
        colnames(relfit) <- c("f=","f>")
        colnames(relcomp) <- c("c=","c>")
        # BF1 <- BF(lm1,hypothesis=paste0("mu_minus_mu0>",as.character(mu0)),
        #           prior=prior)
      } else if(x$alternative=="greater"){
        if(!grepl("Paired",x$method)){
          hypotheses <- c(paste0("true mean is less than or equal to ",as.character(mu0)),paste0("true mean is greater than ",as.character(mu0)))
        }else{
          hypotheses <- c("true difference is smaller than or equal to 0","true difference is greater than 0")
        }
        logBFtu <- c(relfit1-relcomp1,relfit2-relcomp2)
        names(logBFtu) <- hypotheses
        BFtu <- exp(logBFtu)
        PHP <- BFtu/sum(BFtu)
        BFmatrix <- matrix(rep(BFtu,2),ncol=2) / matrix(rep(BFtu,2),ncol=2,byrow=T)
        logBFmatrix <- matrix(rep(logBFtu,2),ncol=2) - matrix(rep(logBFtu,2),ncol=2,byrow=T)
        row.names(BFmatrix) <- row.names(logBFmatrix) <- hypotheses
        colnames(BFmatrix) <- colnames(logBFmatrix) <- hypotheses
        relfit <- matrix(c(rep(1,2),exp(relfit1),exp(relfit2)),ncol=2)
        relcomp <- matrix(c(rep(1,2),rep(.5,2)),ncol=2)
        row.names(relfit) <- row.names(relcomp) <- hypotheses
        colnames(relfit) <- c("f=","f>")
        colnames(relcomp) <- c("c=","c>")
        # BF1 <- BF(lm1,hypothesis=paste0("mu_minus_mu0<",as.character(mu0)),
        #           prior=prior)
      }
    } else if(hypothesis=="exploratory"){
      if(!grepl("Paired",x$method)){
        hypotheses <- c(paste0("true mean is equal to ",as.character(mu0)),paste0("true mean is less than ",as.character(mu0)),paste0("true mean is greater than ",as.character(mu0)))
      }else{
        hypotheses <- c("true difference is equal to 0","true difference is less than 0","true difference is greater than 0")
      }
      logBFtu <- c(relfit0-relcomp0,relfit1-relcomp1,relfit2-relcomp2)
      names(logBFtu) <- hypotheses
      BFtu <- exp(logBFtu)
      PHP <- BFtu/sum(BFtu)
      BFmatrix <- matrix(rep(BFtu,3),ncol=3) / matrix(rep(BFtu,3),ncol=3,byrow=T)
      logBFmatrix <- matrix(rep(logBFtu,3),ncol=3) - matrix(rep(logBFtu,3),ncol=3,byrow=T)
      row.names(BFmatrix) <- row.names(logBFmatrix) <- hypotheses
      colnames(BFmatrix) <- colnames(logBFmatrix) <- hypotheses
      relfit <- matrix(c(exp(relfit0),rep(1,3),exp(relfit1),exp(relfit2)),ncol=2)
      relcomp <- matrix(c(exp(relcomp0),rep(1,3),rep(.5,2)),ncol=2)
      row.names(relfit) <- row.names(relcomp) <- hypotheses
      colnames(relfit) <- c("f=","f>")
      colnames(relcomp) <- c("c=","c>")

      # BF1 <- BF(x=lm1,prior=prior)
    } else {
      stop("hypothesis should be NULL or 'exploratory'.")
    }

  }else{ # two samples t test

    if(!grepl("Welch",x$method)){ # equal variances assumed
      # x=t_test(1:10, y = c(4:8,12)/.4,var.equal = F,alternative="less")
      testname <- x$method
      tvalue <- x$statistic
      mu0 <- x$null.value
      xbar <- x$estimate
      df <- x$parameter
      stderr <- sqrt(sum(x$v*(x$n-1))*sum(1/x$n)/(sum(x$n)-2))
      relfit0 <- dt((xbar[1]-xbar[2])/stderr,df=df,log=TRUE) - log(stderr)
      relfit1 <- log(1-pt((xbar[1]-xbar[2])/stderr,df=df))
      relfit2 <- log(pt((xbar[1]-xbar[2])/stderr,df=df))
      stderr0 <- sqrt(sum((x$v)*(x$n-1)/(x$n)))
      relcomp0 <- dt((xbar[1]-xbar[2])/stderr0,df=1,log=TRUE) - log(stderr0)
      relcomp1 <- relcomp2 <- log(.5)

      if(is.null(hypothesis)){
        if(x$alternative=="two.sided"){
          hypotheses <- c("true difference is equal to 0","true difference is not equal to 0")
          logBFtu <- c(relfit0-relcomp0,0)
          names(logBFtu) <- hypotheses
          BFtu <- exp(logBFtu)
          PHP <- BFtu/sum(BFtu)
          BFmatrix <- matrix(rep(BFtu,2),ncol=2) / matrix(rep(BFtu,2),ncol=2,byrow=T)
          logBFmatrix <- matrix(rep(logBFtu,2),ncol=2) - matrix(rep(logBFtu,2),ncol=2,byrow=T)
          row.names(BFmatrix) <- row.names(logBFmatrix) <- hypotheses
          colnames(BFmatrix) <- colnames(logBFmatrix) <- hypotheses
          relfit <- matrix(c(exp(relfit0),rep(1,3)),ncol=2)
          relcomp <- matrix(c(exp(relcomp0),rep(1,3)),ncol=2)
          row.names(relfit) <- row.names(relcomp) <- hypotheses
          colnames(relfit) <- c("f=","f>")
          colnames(relcomp) <- c("c=","c>")

        }else if(x$alternative=="less"){
          hypotheses <- c("true difference is greater than or equal to 0","true difference is less than 0")
          logBFtu <- c(relfit2-relcomp2,relfit1-relcomp1)
          names(logBFtu) <- hypotheses
          BFtu <- exp(logBFtu)
          PHP <- BFtu/sum(BFtu)
          BFmatrix <- matrix(rep(BFtu,2),ncol=2) / matrix(rep(BFtu,2),ncol=2,byrow=T)
          logBFmatrix <- matrix(rep(logBFtu,2),ncol=2) - matrix(rep(logBFtu,2),ncol=2,byrow=T)
          row.names(BFmatrix) <- row.names(logBFmatrix) <- hypotheses
          colnames(BFmatrix) <- colnames(logBFmatrix) <- hypotheses
          relfit <- matrix(c(rep(1,2),exp(relfit2),exp(relfit1)),ncol=2)
          relcomp <- matrix(c(rep(1,2),rep(.5,2)),ncol=2)
          row.names(relfit) <- row.names(relcomp) <- hypotheses
          colnames(relfit) <- c("f=","f>")
          colnames(relcomp) <- c("c=","c>")
          # BF1 <- BF(lm1,hypothesis=paste0("mu_minus_mu0>",as.character(mu0)),
          #           prior=prior)
        } else if(x$alternative=="greater"){
          hypotheses <- c("true difference is less than or equal to 0","true difference is greater than 0")
          logBFtu <- c(relfit1-relcomp1,relfit2-relcomp2)
          names(logBFtu) <- hypotheses
          BFtu <- exp(logBFtu)
          PHP <- BFtu/sum(BFtu)
          BFmatrix <- matrix(rep(BFtu,2),ncol=2) / matrix(rep(BFtu,2),ncol=2,byrow=T)
          logBFmatrix <- matrix(rep(logBFtu,2),ncol=2) - matrix(rep(logBFtu,2),ncol=2,byrow=T)
          row.names(BFmatrix) <- row.names(logBFmatrix) <- hypotheses
          colnames(BFmatrix) <- colnames(logBFmatrix) <- hypotheses
          relfit <- matrix(c(rep(1,2),exp(relfit1),exp(relfit2)),ncol=2)
          relcomp <- matrix(c(rep(1,2),rep(.5,2)),ncol=2)
          row.names(relfit) <- row.names(relcomp) <- hypotheses
          colnames(relfit) <- c("f=","f>")
          colnames(relcomp) <- c("c=","c>")
          # BF1 <- BF(lm1,hypothesis=paste0("mu_minus_mu0<",as.character(mu0)),
          #           prior=prior)
        }else if(hypothesis=="exploratory"){
          hypotheses <- c("true difference is equal to 0","true difference is less than 0",
                          "true difference is greater than 0")
          logBFtu <- c(relfit0-relcomp0,relfit1-relcomp1,relfit2-relcomp2)
          names(logBFtu) <- hypotheses
          BFtu <- exp(logBFtu)
          PHP <- BFtu/sum(BFtu)
          BFmatrix <- matrix(rep(BFtu,3),ncol=3) / matrix(rep(BFtu,3),ncol=3,byrow=T)
          logBFmatrix <- matrix(rep(logBFtu,3),ncol=3) - matrix(rep(logBFtu,3),ncol=3,byrow=T)
          row.names(BFmatrix) <- row.names(logBFmatrix) <- hypotheses
          colnames(BFmatrix) <- colnames(logBFmatrix) <- hypotheses
          relfit <- matrix(c(exp(relfit0),rep(1,3),exp(relfit1),exp(relfit2)),ncol=2)
          relcomp <- matrix(c(exp(relcomp0),rep(1,3),rep(.5,2)),ncol=2)
          row.names(relfit) <- row.names(relcomp) <- hypotheses
          colnames(relfit) <- c("f=","f>")
          colnames(relcomp) <- c("c=","c>")
        }
      }
    }else{ #equal variances not assumed

      meanN <- x$estimate
      scaleN <- (x$v)/(x$n)
      dfN <- x$n-1
      scale0 <- (x$v)*(x$n-1)/(x$n)
      df0 <- rep(1,2)
      drawsN <- rt(1e5,df=dfN[1])*sqrt(scaleN[1]) + meanN[1] - rt(1e5,df=dfN[2])*sqrt(scaleN[2]) - meanN[2]
      densN <- approxfun(density(drawsN))
      relfit0 <- log(densN(0))
      relfit1 <- log(mean(drawsN<0))
      relfit2 <- log(mean(drawsN>0))
      draws0 <- rt(1e6,df=df0[1])*sqrt(scale0[1])- rt(1e5,df=df0[2])*sqrt(scale0[2])
      relcomp0 <- log((sum((draws0<1)*(draws0> -1))/1e6)/2)
      relcomp1 <- relcomp2 <- log(.5)

      if(is.null(hypothesis)){
        if(x$alternative=="two.sided"){
          hypotheses <- c("true difference is equal to 0","true difference is not equal to 0")
          logBFtu <- c(relfit0-relcomp0,0)
          names(logBFtu) <- hypotheses
          BFtu <- exp(logBFtu)
          PHP <- BFtu/sum(BFtu)
          BFmatrix <- matrix(rep(BFtu,2),ncol=2) / matrix(rep(BFtu,2),ncol=2,byrow=T)
          logBFmatrix <- matrix(rep(logBFtu,2),ncol=2) - matrix(rep(logBFtu,2),ncol=2,byrow=T)
          row.names(BFmatrix) <- row.names(logBFmatrix) <- hypotheses
          colnames(BFmatrix) <- colnames(logBFmatrix) <- hypotheses
          relfit <- matrix(c(exp(relfit0),rep(1,3)),ncol=2)
          relcomp <- matrix(c(exp(relcomp0),rep(1,3)),ncol=2)
          row.names(relfit) <- row.names(relcomp) <- hypotheses
          colnames(relfit) <- c("f=","f>")
          colnames(relcomp) <- c("c=","c>")

        }else if(x$alternative=="less"){
          hypotheses <- c("true difference is greater than or equal to 0","true difference is less than 0")
          logBFtu <- c(relfit2-relcomp2,relfit1-relcomp1)
          names(logBFtu) <- hypotheses
          BFtu <- exp(logBFtu)
          PHP <- BFtu/sum(BFtu)
          BFmatrix <- matrix(rep(BFtu,2),ncol=2) / matrix(rep(BFtu,2),ncol=2,byrow=T)
          logBFmatrix <- matrix(rep(logBFtu,2),ncol=2) - matrix(rep(logBFtu,2),ncol=2,byrow=T)
          row.names(BFmatrix) <- row.names(logBFmatrix) <- hypotheses
          colnames(BFmatrix) <- colnames(logBFmatrix) <- hypotheses
          relfit <- matrix(c(rep(1,2),exp(relfit2),exp(relfit1)),ncol=2)
          relcomp <- matrix(c(rep(1,2),rep(.5,2)),ncol=2)
          row.names(relfit) <- row.names(relcomp) <- hypotheses
          colnames(relfit) <- c("f=","f>")
          colnames(relcomp) <- c("c=","c>")
          # BF1 <- BF(lm1,hypothesis=paste0("mu_minus_mu0>",as.character(mu0)),
          #           prior=prior)
        } else if(x$alternative=="greater"){
          hypotheses <- c("true difference is less than or equal to 0","true difference is greater than 0")
          logBFtu <- c(relfit1-relcomp1,relfit2-relcomp2)
          names(logBFtu) <- hypotheses
          BFtu <- exp(logBFtu)
          PHP <- BFtu/sum(BFtu)
          BFmatrix <- matrix(rep(BFtu,2),ncol=2) / matrix(rep(BFtu,2),ncol=2,byrow=T)
          logBFmatrix <- matrix(rep(logBFtu,2),ncol=2) - matrix(rep(logBFtu,2),ncol=2,byrow=T)
          row.names(BFmatrix) <- row.names(logBFmatrix) <- hypotheses
          colnames(BFmatrix) <- colnames(logBFmatrix) <- hypotheses
          relfit <- matrix(c(rep(1,2),exp(relfit1),exp(relfit2)),ncol=2)
          relcomp <- matrix(c(rep(1,2),rep(.5,2)),ncol=2)
          row.names(relfit) <- row.names(relcomp) <- hypotheses
          colnames(relfit) <- c("f=","f>")
          colnames(relcomp) <- c("c=","c>")
          # BF1 <- BF(lm1,hypothesis=paste0("mu_minus_mu0<",as.character(mu0)),
          #           prior=prior)
        }else if(hypothesis=="exploratory"){
          hypotheses <- c("true difference is equal to 0","true difference is less than 0",
                          "true difference is greater than 0")
          logBFtu <- c(relfit0-relcomp0,relfit1-relcomp1,relfit2-relcomp2)
          names(logBFtu) <- hypotheses
          BFtu <- exp(logBFtu)
          PHP <- BFtu/sum(BFtu)
          BFmatrix <- matrix(rep(BFtu,3),ncol=3) / matrix(rep(BFtu,3),ncol=3,byrow=T)
          logBFmatrix <- matrix(rep(logBFtu,3),ncol=3) - matrix(rep(logBFtu,3),ncol=3,byrow=T)
          row.names(BFmatrix) <- row.names(logBFmatrix) <- hypotheses
          colnames(BFmatrix) <- colnames(logBFmatrix) <- hypotheses
          relfit <- matrix(c(exp(relfit0),rep(1,3),exp(relfit1),exp(relfit2)),ncol=2)
          relcomp <- matrix(c(exp(relcomp0),rep(1,3),rep(.5,2)),ncol=2)
          row.names(relfit) <- row.names(relcomp) <- hypotheses
          colnames(relfit) <- c("f=","f>")
          colnames(relcomp) <- c("c=","c>")
        }
      }
    }
  }
  return(list(BFtu=BFtu,logBFtu=logBFtu,PHP=PHP,BFmatrix=BFmatrix,logBFmatrix=logBFmatrix,
              relfit=relfit,relcomp=relcomp))
}


BFupdate.htest <- function(BF1,
                     x,
                     ...){
  stop("REMINDER: Needs to be done.")
}




