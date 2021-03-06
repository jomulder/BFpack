

#' @method BF hetcor
#' @export
BF.hetcor <- function(x,
                       hypothesis = NULL,
                       prior.hyp = NULL,
                       ...){
  get_est <- get_estimates(x)
  P <- nrow(x$std.errors)
  numcorr <- P*(P-1)/2
  estimates <- get_est$estimate
  errcov <- get_est$Sigma[[1]]

  corr_names <- names(get_est$estimate)
  matrix_names <- matrix(corr_names,nrow=P)

  #exploratory BF testing
  relfit <- matrix(c(dnorm(0,mean=estimates,sd=sqrt(diag(errcov))),
                     pnorm(0,mean=estimates,sd=sqrt(diag(errcov))),
                     1-pnorm(0,mean=estimates,sd=sqrt(diag(errcov)))),ncol=3)
  relcomp <- matrix(c(dbeta(rep(.5,numcorr),shape1=P/2,shape2=P/2)*.5,
                      rep(.5,2*numcorr)),ncol=3)
  row.names(relfit) <- row.names(relcomp) <- names(estimates)

  BFtu_exploratory <- relfit / relcomp
  colnames(BFtu_exploratory) <- colnames(BFtu_exploratory) <-  c("Pr(=0)","Pr(<0)","Pr(>0)")
  PHP_exploratory <- BFtu_exploratory / apply(BFtu_exploratory,1,sum)

  #confirmatory BF testing
  if(!is.null(hypothesis)){
    parse_hyp <- parse_hypothesis(names(estimates),hypothesis)
    parse_hyp$hyp_mat <- do.call(rbind, parse_hyp$hyp_mat)
    #create coefficient with equality and order constraints
    RrList <- make_RrList2(parse_hyp)
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]

    numhyp <- length(RrE)
    relcomp <- t(matrix(unlist(
      lapply(1:numhyp, function(h){
        jointuniform_measures(P,numcorr, 1, RrE[[h]], RrO[[h]], Fisher=0)
    })
    ),nrow=2))
    relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Gaussian_measures(estimates,errcov,RrE1=RrE[[h]],RrO1=RrO[[h]],names1=names(estimates),
                        constraints1=parse_hyp$original_hypothesis[h])
    })),nrow=2))
    row.names(relcomp) <- parse_hyp$original_hypothesis
    row.names(relfit) <- parse_hyp$original_hypothesis
    # evaluation of complement hypothesis
    relfit <- Gaussian_prob_Hc(estimates,errcov,relfit,RrO)
    relcomp <- jointuniform_prob_Hc(P,numcorr,1,relcomp,RrO)

    colnames(relcomp) <- c("c_E","c_O")
    colnames(relfit) <- c("f_E","f_O")
    # computation of exploratory BFs and PHPs
    # the BF for the complement hypothesis vs Hu needs to be computed.
    BFtu_confirmatory <- c(apply(relfit / relcomp, 1, prod))
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
    PHP_confirmatory <- BFtu_confirmatory*priorprobs / sum(BFtu_confirmatory*priorprobs)
    BFtable <- cbind(relcomp,relfit,relfit[,1]/relcomp[,1],relfit[,2]/relcomp[,2],
                     apply(relfit,1,prod)/apply(relcomp,1,prod),PHP_confirmatory)
    row.names(BFtable) <- names(BFtu_confirmatory)
    colnames(BFtable) <- c("complex=","complex>","fit=","fit>","BF=","BF>","BF","PHP")
    BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory))/
      t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
    diag(BFmatrix_confirmatory) <- 1
    row.names(BFmatrix_confirmatory) <- colnames(BFmatrix_confirmatory) <- names(BFtu_confirmatory)
    hypotheses <- row.names(relcomp)
  }else{
    BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relfit <-
      relcomp <- hypotheses <- BFtable <- priorprobs <- NULL
  }

  # Store in output
  BF_out <- list(
    BFtu_exploratory=BFtu_exploratory,
    PHP_exploratory=PHP_exploratory,
    BFtu_confirmatory=BFtu_confirmatory,
    PHP_confirmatory=PHP_confirmatory,
    BFmatrix_confirmatory=BFmatrix_confirmatory,
    BFtable_confirmatory=BFtable,
    prior.hyp=priorprobs,
    hypotheses=hypotheses,
    estimates=estimates,
    model=x,
    bayesfactor="Bayes factor based on joint uniform priors",
    parameter="Correlations",
    call=match.call())

  class(BF_out) <- "BF"

  return(BF_out)

}

