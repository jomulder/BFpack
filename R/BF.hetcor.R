

#' @method BF hetcor
#' @export
BF.hetcor <- function(x,
                       hypothesis = NULL,
                       prior.hyp.explo = NULL,
                       prior.hyp.conf = NULL,
                       prior.hyp = NULL,
                       complement = TRUE,
                       log = FALSE,
                       cov.prob = .95,
                       ...){

  if(!(cov.prob>0 & cov.prob<1)){
    stop("The argument 'cov.prob' is a coverage probability for the interval estimates that should lie between 0 and 1. The default is 0.95.")
  }
  CrI_LB <- (1 - cov.prob)/2
  CrI_UB <- 1 - (1 - cov.prob)/2

  logIN <- log

  # check proper usage of argument 'prior.hyp.conf' and 'prior.hyp.explo'
  if(!is.null(prior.hyp.conf)){
    prior.hyp <- prior.hyp.conf
  }
  prior.hyp.explo <- process.prior.hyp.explo(prior_hyp_explo = prior.hyp.explo, model=x)

  get_est <- get_estimates(x)
  P <- nrow(x$std.errors)
  numcorr <- P*(P-1)/2
  estimates <- get_est$estimate
  errcov <- get_est$Sigma[[1]]
  # use Fisher transformed for both exploratory and confirmatory to get consistent results.
  # skewness in the likelihood is ignored.
  est.var.F <- do.call(cbind,lapply(1:numcorr,function(c){
    draws.norm <- rnorm(1e5,mean=estimates[c],sd=sqrt(errcov[c,c]))
    draws.norm.F <- FisherZ(draws.norm[draws.norm < 1 & draws.norm > -1])
    return(c(median(draws.norm.F),var(draws.norm.F)))
  }))
  estimates.F <- est.var.F[1,]
  if(numcorr > 1){
    errcov.F <- diag(est.var.F[2,])
  }else{
    errcov.F <- as.matrix(est.var.F[2,])
  }
  corr_names <- names(get_est$estimate)
  #matrix_names <- matrix(corr_names,nrow=P)
  names(estimates.F) <- colnames(errcov.F) <- row.names(errcov.F) <- corr_names

  #exploratory BF testing
  relfit <- matrix(c(dnorm(0,mean=estimates.F,sd=sqrt(diag(errcov.F)),log=TRUE),
                     pnorm(0,mean=estimates.F,sd=sqrt(diag(errcov.F)),log.p=TRUE),
                     pnorm(0,mean=estimates.F,sd=sqrt(diag(errcov.F)),log.p=TRUE,
                           lower.tail=FALSE)),ncol=3)
  # get draws from joint uniform prior to compute relative measures
  if(sum(P==Fcor$P)==0){
    #number of draws to get 1e7 draws for the marginal of 1 Fisher transformation correlation
    numdraws <- round(1e7/(P*(P-1)/2))
    drawsJU <- draw_ju_r(P,samsize=numdraws,Fisher=1)
    approx_studt <- QRM::fit.st(c(drawsJU))$par.ests[c(1,3)]
  }else{
    # use the estimate of the scale from the Fcor object
    # for the df the estimates show some numerical error for P>20, so use fitted line
    approx_studt <- unlist(c(Fcor[which(P==Fcor$P),1:2]))
    if(P > 20){
      # use fitted linear line rather than rough estimate of df
      slpe1 <- 2.944494
      intcept1 <- 1.864901
      approx_studt[1] <- P * slpe1 + intcept1
    }
  }
  relcomp0 <- dt(0,df=approx_studt[1],log=TRUE)-log(approx_studt[2]) # all marginal priors are the same
  relcomp <- matrix(c(rep(relcomp0,numcorr),rep(log(.5),numcorr*2)),ncol=3)
  row.names(relfit) <- row.names(relcomp) <- names(estimates.F)

  BFtu_exploratory <- relfit - relcomp
  colnames(BFtu_exploratory) <- colnames(BFtu_exploratory) <-  c("Pr(=0)","Pr(<0)","Pr(>0)")
  maxrow <- apply(BFtu_exploratory,1,max)
  BFtu_explo_norm <- exp(BFtu_exploratory - maxrow %*% t(rep(1,ncol(BFtu_exploratory)))) *
    rep(1,nrow(BFtu_exploratory)) %*% t(prior.hyp.explo[[1]])
  PHP_exploratory <- BFtu_explo_norm / apply(BFtu_explo_norm,1,sum)

  if(logIN == FALSE){
    BFtu_exploratory <- exp(BFtu_exploratory)
  }

  #confirmatory BF testing
  if(!is.null(hypothesis)){
    numG <- 1
    numcorrgroup <- numcorr
    varnames <- list(row.names(x$correlations))
    # get all names combinations for correlations (similar as BF.cor_test)
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
    x$corrnames <- corrnames

    params_in_hyp1 <- params_in_hyp(hypothesis)

    corr_names <- unlist(lapply(1:length(x$corrnames),function(g){
      c(x$corrnames[[g]][lower.tri(x$corrnames[[g]])],
        t(x$corrnames[[g]])[lower.tri(x$corrnames[[g]])])
    })) #which includes Y1_with_Y2 and Y2_with_Y1

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

    numhyp <- length(RrE)
    relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Gaussian_measures(estimates,errcov,RrE1=RrE[[h]],RrO1=RrO[[h]],names1=names(estimates),
                        constraints1=parse_hyp$original_hypothesis[h])
    })),nrow=2))
    # approximate unconstrained Fisher transformed correlations with a multivariate Student t
    mean0 <- rep(0,numcorr)
    if(numcorr==1){
      Scale0 <- as.matrix(approx_studt[2]**2)
      df0 <- round(approx_studt[1])
    }else{
      Scale0 <- diag(rep(approx_studt[2]**2,numcorr))
      df0 <- round(approx_studt[1])
    }
    mean0 <- rep(0,numcorr)
    relcomp <- t(matrix(unlist(lapply(1:numhyp,function(h){
      relcomp_h <- Student_measures(mean1=mean0,
                                    Scale1=Scale0,
                                    df1=df0,
                                    RrE1=RrE[[h]],
                                    RrO1=RrO[[h]])
      return(relcomp_h)

    })),nrow=2))

    row.names(relcomp) <- parse_hyp$original_hypothesis
    row.names(relfit) <- parse_hyp$original_hypothesis
    # evaluation of complement hypothesis
    if(complement == TRUE){
      relfit <- Gaussian_prob_Hc(estimates.F,errcov.F,relfit,RrO)
      relcomp <- Student_prob_Hc(mean1=mean0,scale1=Scale0,df1=df0,relmeas1=relcomp,
                                 constraints=NULL,RrO1=RrO)
    }
    hypothesisshort <- unlist(lapply(1:nrow(relfit),function(h) paste0("H",as.character(h))))
    row.names(relfit) <- row.names(relfit) <- hypothesisshort

    colnames(relcomp) <- c("c_E","c_O")
    colnames(relfit) <- c("f_E","f_O")
    # computation of exploratory BFs and PHPs
    # the BF for the complement hypothesis vs Hu needs to be computed.
    BFtu_confirmatory <- c(apply(relfit - relcomp, 1, sum))
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
    maxBFtu <- max(BFtu_confirmatory)
    PHP_confirmatory <- exp(BFtu_confirmatory-maxBFtu)*priorprobs /
      sum(exp(BFtu_confirmatory-maxBFtu)*priorprobs)
    BFtable <- cbind(relcomp,relfit,relfit[,1]-relcomp[,1],relfit[,2]-relcomp[,2],
                     apply(relfit,1,sum)-apply(relcomp,1,sum),PHP_confirmatory)
    BFtable[,1:7] <- exp(BFtable[,1:7])
    row.names(BFtable) <- names(BFtu_confirmatory)
    colnames(BFtable) <- c("complex=","complex>","fit=","fit>","BF=","BF>","BF","PHP")
    BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)) -
      t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
    diag(BFmatrix_confirmatory) <- 0
    row.names(BFmatrix_confirmatory) <- colnames(BFmatrix_confirmatory) <- names(BFtu_confirmatory)
    hypotheses <- row.names(relcomp)

    if(logIN == FALSE){
      BFtu_confirmatory <- exp(BFtu_confirmatory)
      BFmatrix_confirmatory <- exp(BFmatrix_confirmatory)
    }

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
    bayesfactor="Bayes factors based on joint uniform priors",
    parameter="measures of association",
    log=logIN,
    call=match.call())

  class(BF_out) <- "BF"

  return(BF_out)

}
