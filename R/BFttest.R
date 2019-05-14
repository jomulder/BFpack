

BF.htest <- function(x,
                      hypothesis = NULL,
                      prior = NULL,
                      ...){

  numpop <- length(x$estimate)

  if(numpop==1){ #one sample t test
    tvalue <- x$statistic
    names(tvalue) <- NULL
    n <- x$parameter + 1
    names(n) <- NULL
    mu0 <- x$null.value
    names(mu0) <- NULL
    xbar <- x$estimate
    names(xbar) <- NULL
    sd1 <- sqrt(n) * (xbar - mu0) / tvalue
    y1 <- rnorm(n)
    y1 <- y1 - mean(y1)
    y1 <- sd1*y1/sd(y1) + xbar - mu0
    lm1 <- lm(y1 ~ 1)
    names(lm1$coefficients) <- "mu_minus_mu0"

    if(is.null(hypothesis)){
      if(x$alternative=="two.sided"){
        BF1 <- BF(lm1,hypothesis=paste0("mu_minus_mu0=",as.character(mu0)),
                  prior=prior)
      } else if(x$alternative=="less"){
        BF1 <- BF(lm1,hypothesis=paste0("mu_minus_mu0<",as.character(mu0)),
                  prior=prior)
      } else if(x$alternative=="greater"){
        BF1 <- BF(lm1,hypothesis=paste0("mu_minus_mu0>",as.character(mu0)),
                  prior=prior)
      }
    } else if(hypothesis=="exploratory"){
      BF1 <- BF(x=lm1,prior=prior)
    } else {
      if(!grepl("mu_minus_mu0",hypothesis)){stop("parameter in 'hypothesis' to be tested calls 'mu_minus_mu0'.")}
      BF1 <- BF(lm1,hypothesis=hypothesis,prior=prior)
    }

  }else{ # two samples t test
    stop("for a two sample t test please use an aov or lm object instead of a htest object.")
  }

  return(BF1)

}


BFupdate.htest <- function(BF1,
                     x,
                     ...){
  stop("REMINDER: Needs to be done.")
}




