

#' @importFrom stats terms
remove_predictors_helper <- function(Y_groups, formula){

  # number of groups
  groups <- length(Y_groups)

  # model matrix terms
  mm_terms <- attr(terms(formula), "term.labels")

  if(length(mm_terms) == 0){

    Y_groups

  } else {
    lapply(seq_len(groups), function(x){

      # check for factors
      factor_pred <- which(paste0("as.factor(", colnames(Y_groups[[x]]), ")") %in% mm_terms)

      # check for non factors
      cont_pred <- which(colnames(Y_groups[[x]]) %in% mm_terms)

      # remove predictors
      Y_groups[[x]][,-c(factor_pred, cont_pred)]

    })
  }
}


FisherZ <- function(r){.5*log((1+r)/(1-r))}


globalVariables(c("Fcor"))

#' @importFrom mvtnorm dmvnorm pmvnorm rmvnorm
#' @importFrom utils globalVariables
#' @importFrom stats dnorm pnorm
#' @importFrom QRM fit.st
#' @method BF cor_test
#' @export
BF.cor_test <- function(x,
                        hypothesis = NULL,
                        prior.hyp.explo = NULL,
                        prior.hyp.conf = NULL,
                        prior.hyp = NULL,
                        complement = TRUE,
                        log = FALSE,
                        cov.prob = .95,
                        ...){

  if(x$prior.cor == "joint.unif"){
    bayesfactor <- "Bayes factors based on joint uniform priors"
  }
  if(x$prior.cor == "marg.unif"){
    bayesfactor <- "Bayes factors based on marginally uniform priors"
  }
  if( !(x$prior.cor == "marg.unif" | x$prior.cor == "joint.unif") ){
    stop("argument 'prior' of cor_test object should be 'joint.unif' or 'marg.unif'. See ?cor_test")
  }
  testedparameter <- "correlation coefficients"

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
  prior.hyp.explo <- process.prior.hyp.explo(prior_hyp_explo = prior.hyp.explo, model = x)

  P <- dim(x$corrdraws[[1]])[2]
  numG <- length(x$corrdraws)
  numcorrgroup <- P*(P-1)/2
  get_est <- get_estimates(x)
  corrmeanN <- get_est$estimate
  corrcovmN <- get_est$Sigma[[1]]

  # Exploratory testing of correlation coefficients
  # get height of prior density at 0 of Fisher transformed correlation
  # slope and intercept of df as function of P based on Fcor

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
  if(x$prior.cor == "joint.unif"){
    relcomp0 <- dt(0,df=approx_studt[1],log=TRUE)-log(approx_studt[2]) # all marginal priors are the same
  }
  if(x$prior.cor == "marg.unif"){
    relcomp0 <- dt(0,df=Fcor[1,1],log=TRUE)-log(Fcor[1,2]) #Fisher transformed height for a uniform(-1,1) prior
  }

  # compute exploratory BFs
  corr_names <- rownames(x$correstimates)
  numcorr <- length(corrmeanN)
  relfit <- matrix(c(dnorm(0,mean=corrmeanN,sd=sqrt(diag(corrcovmN)),log=TRUE),
                     pnorm(0,mean=corrmeanN,sd=sqrt(diag(corrcovmN)),log.p=TRUE),
                     pnorm(0,mean=corrmeanN,sd=sqrt(diag(corrcovmN)),log.p=TRUE,
                             lower.tail=FALSE)),ncol=3)
  relcomp <- matrix(c(rep(relcomp0,numcorr),rep(log(.5),numcorr*2)),ncol=3)
  colnames(relcomp) <- colnames(relfit) <- c("p(=0)","Pr(<0)","Pr(>0)")
  BFtu_exploratory <- relfit - relcomp
  row.names(BFtu_exploratory) <- rownames(x$correstimates)
  colnames(BFtu_exploratory) <- c("BF0u","BF1u","BF2u")
  rowmax <- apply(BFtu_exploratory,1,max)
  norm_BF_explo <- exp(BFtu_exploratory - rowmax %*% t(rep(1,ncol(BFtu_exploratory)))) *
    (rep(1,nrow(BFtu_exploratory)) %*% t(prior.hyp.explo[[1]]))
  PHP_exploratory <- norm_BF_explo / apply(norm_BF_explo,1,sum)
  colnames(PHP_exploratory) <- c("P(=0)","P(<0)","P(>0)")

  # posterior estimates
  postestimates_correlations <- Reduce(rbind,
                                       lapply(1:numG,function(g){
                                         draws_stack_g <- as.matrix(do.call(cbind,lapply(1:P,function(p){x$corrdraws[[g]][,,p]}))[,which(c(lower.tri(diag(P))))])
                                         means <- apply(draws_stack_g,2,mean)
                                         medians <- apply(draws_stack_g,2,median)
                                         lb <- apply(draws_stack_g,2,quantile,CrI_LB)
                                         ub <- apply(draws_stack_g,2,quantile,CrI_UB)
                                         rm(draws_stack_g)
                                         return(cbind(means,medians,lb,ub))
                                       }))

  colnames(postestimates_correlations) <- c("mean","median",paste0(as.character(round(CrI_LB*100,7)),"%"),
                                            paste0(as.character(round(CrI_UB*100,7)),"%"))
  postestimates <- postestimates_correlations
  rownames(postestimates) <- corr_names

  if(logIN == FALSE){
    BFtu_exploratory <- exp(BFtu_exploratory)
  }

  # confirmatory testing if hypothesis argument is used
  if(!is.null(hypothesis)){
    if(x$prior.cor != "joint.unif"){
      BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relfit <-
        relcomp <- hypotheses <- BFtable <- priorprobs <- NULL
      warning(
"'hypothesis' argument not used. No manual hypothesis test is supported when the prior of the
correlation matrix is marginally uniform. Set the 'prior.cor' argument in cor_test() to 'joint.unif'
for manual hypothesis testing (Mulder and Gelissen, 2023).")
    }else{

      #check if constraints are formulated on correlations in different populations
      #if so, then the correlation names contains the string "_in_g" at the end
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
      #Fisher transform constants in constraints
      for(h in 1:numhyp){
        if(!is.null(RrE[[h]])){
          RrE[[h]][,ncol(RrE[[h]])] <- FisherZ(RrE[[h]][,ncol(RrE[[h]])])
        }
        if(!is.null(RrO[[h]])){
          RrO[[h]][,ncol(RrO[[h]])] <- FisherZ(RrO[[h]][,ncol(RrO[[h]])])
        }
      }

      relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
        Gaussian_measures(corrmeanN,corrcovmN,RrE1=RrE[[h]],RrO1=RrO[[h]])
      })),nrow=2))
      #names1 and constraints1 ... to fix ...
      # approximate unconstrained Fisher transformed correlations with a multivariate Student t
      if(numcorrgroup==1){
        if(numcorr==1){
          Scale0 <- as.matrix(approx_studt[2]**2)
        }else{
          Scale0 <- diag(rep(approx_studt[2]**2,numG))
        }
        mean0 <- rep(0,numG)
        df0 <- round(approx_studt[1])
      }else{
        mean0 <- rep(0,numcorrgroup*numG)
        Scale0 <- diag(rep(approx_studt[2]**2,numcorrgroup*numG))
        df0 <- round(approx_studt[1])
      }
      relcomp <- t(matrix(unlist(lapply(1:numhyp,function(h){
        relcomp_h <- Student_measures(mean1=mean0,
                                      Scale1=Scale0,
                                      df1=df0,
                                      RrE1=RrE[[h]],
                                      RrO1=RrO[[h]])
        return(relcomp_h)

      })),nrow=2))

      row.names(relfit) <- row.names(relcomp) <- parse_hyp$original_hypothesis
      if(complement == TRUE){
        relfit <- Gaussian_prob_Hc(corrmeanN,corrcovmN,relfit,RrO)
        relcomp <- Student_prob_Hc(mean1=mean0,scale1=Scale0,df1=df0,relmeas1=relcomp,
                                   constraints=NULL,RrO1=RrO)
      }
      hypothesisshort <- unlist(lapply(1:nrow(relfit),function(h) paste0("H",as.character(h))))
      row.names(relfit) <- row.names(relfit) <- hypothesisshort

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
      PHP_confirmatory <- exp(BFtu_confirmatory-maxBFtu)*priorprobs / sum(exp(BFtu_confirmatory-maxBFtu)*priorprobs)
      BFtable <- cbind(relcomp,relfit,relfit[,1]-relcomp[,1],relfit[,2]-relcomp[,2],
                       apply(relfit,1,sum)-apply(relcomp,1,sum),PHP_confirmatory)
      BFtable[,1:7] <- exp(BFtable[,1:7])
      row.names(BFtable) <- names(BFtu_confirmatory)
      colnames(BFtable) <- c("complex=","complex>","fit=","fit>","BF=","BF>","BF","PHP")
      BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory))-
        t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
      diag(BFmatrix_confirmatory) <- 0
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
    }
  }else{
    BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relfit <-
      relcomp <- hypotheses <- BFtable <- priorprobs <- NULL
  }

  BFcorr_out <- list(
    BFtu_exploratory=BFtu_exploratory,
    PHP_exploratory=PHP_exploratory,
    BFtu_confirmatory=BFtu_confirmatory,
    PHP_confirmatory=PHP_confirmatory,
    BFmatrix_confirmatory=BFmatrix_confirmatory,
    BFtable_confirmatory=BFtable,
    prior.hyp.explo=prior.hyp.explo,
    prior.hyp.conf=priorprobs,
    hypotheses=hypotheses,
    estimates=postestimates,
    model=x,
    bayesfactor=bayesfactor,
    parameter=testedparameter,
    log=logIN,
    call=match.call()
  )

  class(BFcorr_out) <- "BF"

  return(BFcorr_out)

}


#get draws from joint uniform prior in Fisher transformed space
#Call Fortran subroutine in from bct_prior.f90
draw_ju_r <- function(P, samsize=50000, Fisher=1){
  testm <- matrix(0,ncol=.5*P*(P-1),nrow=samsize)
  #  random1 <- rnorm(1)
  #  random1 <- (random1 - floor(random1))*1e6
  res <-.Fortran("draw_ju",P = as.integer(P),
                 drawscorr=testm,
                 samsize=as.integer(samsize),
                 numcorrgroup=as.integer(.5*P*(P-1)),
                 Fisher=as.integer(Fisher),
                 #seed=as.integer( sample.int(1e6,1) ),
                 PACKAGE="BFpack")
  return(res$drawscorr)

}


#' @title Bayesian correlation analysis
#'
#' @name cor_test
#'
#' @description Estimate the unconstrained posterior for the correlations using a joint uniform prior (Mulder and Gelissen,
#' 2023) or a marginally uniform prior (Barnard et al., 2000, Mulder, 2016).
#'
#' @param ...  matrices (or data frames) of dimensions \emph{n} (observations) by  \emph{p} (variables)
#' for different groups (in case of multiple matrices or data frames).
#'
#' @param formula an object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (e.g., \code{~ education}).
#'
#' @param iter number of iterations from posterior (default is 5000).
#'
#' @param burnin number of iterations for burnin (default is 3000).
#'
#' @param nugget.scale a scalar to avoid computational issues due to posterior draws for the corralations
#' too close to 1 in absolute value. Posterior draws for the correlations are multiplied with this nugget.scale.
#' So \code{nugget.scale} should be close to 1 (the default is .995). If the traceplots show that draws are stuck
#' at 1 or -1 too long try a slightly smaller \code{nugget.scale}.
#'
#' @param prior.cor setting this argument to \code{joint.unif} uses the joint uniform prior for the correlation matrix
#' (Mulder and Gelissen, 2023) and setting this to \code{marg.unif} implies a marginal uniform prior (Barnard et al., 2000,
#' Mulder, 2016).
#'
#' @return list of class \code{cor_test}:
#' \itemize{
#' \item \code{meanF} posterior means of Fisher transform correlations
#' \item \code{covmF} posterior covariance matrix of Fisher transformed correlations
#' \item \code{correstimates} posterior estimates of correlation coefficients
#' \item \code{corrdraws} list of posterior draws of correlation matrices per group
#' \item \code{corrnames} names of all correlations
#' }
#'
#' @references Mulder, J., & Gelissen, J. P. (2023). Bayes factor testing of equality and order constraints on measures of
#' association in social research. Journal of Applied Statistics, 50(2), 315-351. <https://doi.org/10.1080/02664763.2021.1992360>
#'
#' @references Mulder, J. (2016). Bayes factors for testing order-constrained hypotheses on correlations. Journal of Mathematical
#' Psychology, 72, 104-115. <https://doi.org/10.1016/j.jmp.2014.09.004>
#'
#' @references Barnard, J., McCulloch, R., & Meng, X. L. (2000). Modeling covariance matrices in terms of standard deviations and
#' correlations, with application to shrinkage. Statistica Sinica, 1281-1311. <https://www.jstor.org/stable/24306780>
#'
#' @references Joe. Generating random correlation matrices based on partial correlations, Journal of Multivariate Analysis,
#' 97(10), 2177-2189. <https://doi.org/10.1016/j.jmva.2005.05.010>
#'
#' @examples
#' \donttest{
#' # Bayesian correlation analysis of the 6 variables in 'memory' object
#' # we consider a correlation analysis of the first three variable of the memory data.
#' fit <- cor_test(BFpack::memory[,1:3])
#'
#' # Bayesian correlation of variables in memory object in BFpack while controlling
#' # for the Cat variable
#' fit <- cor_test(BFpack::memory[,c(1:4)],formula = ~ Cat)
#'
#' # Example of Bayesian estimation of polyserial correlations
#' memory_example <- memory[,c("Im","Rat")]
#' memory_example$Rat <- as.ordered(memory_example$Rat)
#' fit <- cor_test(memory_example)
#'
#' # Bayesian correlation analysis of first three variables in memory data
#' # for two different groups
#' HC <- subset(BFpack::memory[,c(1:3,7)], Group == "HC")[,-4]
#' SZ <- subset(BFpack::memory[,c(1:3,7)], Group == "SZ")[,-4]
#' fit <- cor_test(HC,SZ)
#'
#' }
#' @rdname cor_test
#' @export
cor_test <- function(..., formula = NULL, iter = 5e3, burnin = 3e3, nugget.scale = .995, prior.cor = "joint.unif"){

  if(is.na(prior.cor)){stop("'prior.cor' argument needs to be either 'joint.unif' or 'marg.unif'. See ?cor_test.")}
  if(is.null(prior.cor)){stop("'prior.cor' argument needs to be either 'joint.unif' or 'marg.unif'. See ?cor_test.")}
  if(!(prior.cor == "joint.unif" | prior.cor == "marg.unif")){
    stop("'prior.cor' argument needs to be either 'joint.unif' or 'marg.unif'. See ?cor_test.")}
  if(prior.cor == "joint.unif"){
    priorchoice <- 1
  }else{
    priorchoice <- 2
  }

  if(!is.numeric(nugget.scale)){stop("'nugget.scale' should be a numerical scalar.")}
  if(nugget.scale > 1 | nugget.scale < 0){stop("'nugget.scale' should be very close 1. If should not exceed 1 nor fall below 0.")}
  nugget.scale <- nugget.scale[1]

  Y_input <- list(...)
  numG <- length(Y_input)

  if(is.null(formula)){
    formula <- ~ 1
  }
  Xnames <- attr(terms(formula), "term.labels")
  whichDV <- lapply(Y_input,function(y){
    unlist(lapply(colnames(y),function(x){sum(x==Xnames)==0}))
  })
  if(numG>1){ #check that the same number of DVs are present in each group (that's how dimensions are coded)
    numDV <- rep(NA,numG)
    for(gg in 1:numG){
      numDV[gg] <- sum(whichDV[[gg]])
    }
    if(sum(abs(diff(numDV)))!=0){
      stop("Each group should contain same number of dependent variables.")
    }
  }

  #check measurement level of dependent variables
  P <- sum(whichDV[[1]])
  ordi <- numcats <- matrix(0,nrow=numG,ncol=P)
  Ylevel <- matrix(0,nrow=numG,ncol=P)
  Y_groups <- Y_input
  for(gg in 1:numG){
    teller <- 1
    for(pp in which(whichDV[[gg]])){
      if(class(Y_groups[[gg]][,pp])[1] == "numeric" | class(Y_groups[[gg]][,pp])[1] == "integer"){
        teller <- teller + 1
        Ylevel[gg,pp] <- "numeric"
      }else{
        if(class(Y_groups[[gg]][,pp])[1] == "ordered"){
          #levels(Y_groups[[gg]][,pp]) <- 1:length(levels(Y_groups[[gg]][,pp]))
          old_levels <- sort(as.numeric(unique(Y_groups[[gg]][,pp])))
          for(index_levels_g_p in 1:length(old_levels)){
            Y_groups[[gg]][which(Y_groups[[gg]][,pp] == old_levels[index_levels_g_p]),pp] <- index_levels_g_p
          }
          Y_groups[[gg]][,pp] <- as.numeric(Y_groups[[gg]][,pp])
          ordi[gg,teller] <- 1
          numcats[gg,teller] <- max(Y_groups[[gg]][,pp])
          Ylevel[gg,pp] <- "ordinal"
          if(numcats[gg,teller]==2){Ylevel[gg,pp] <- "dichotomous"}
          teller <- teller + 1
          if(max(Y_groups[[gg]][,pp])>11){
            stop("Ordinal variables are not allowed to have more than 11 categories")
          }
        }else{
          if(class(Y_groups[[gg]][,pp])[1] == "factor"){
            if(length(levels(Y_groups[[gg]][,pp]))==2){
              Ylevel[gg,pp] <- "dichotomous"
              levels(Y_groups[[gg]][,pp]) <- 1:length(levels(Y_groups[[gg]][,pp]))
              Y_groups[[gg]][,pp] <- as.numeric(Y_groups[[gg]][,pp])
              ordi[gg,teller] <- 1
              numcats[gg,teller] <- 2
              teller <- teller + 1
            }else{
              stop("Outcome variables should be either of class 'numeric', 'ordered', or a 2-level 'factor'.")
            }
          }else{
            stop("Outcome variables should be either of class 'numeric', 'ordered', or a 2-level 'factor'.")
          }
        }
      }
    }
    ordi[gg,] <- as.double(ordi[gg,])
    numcats[gg,] <- as.double(numcats[gg,])
  }
  if(max(numcats) == 1){
    stop("One categorical variable is constant.")
  }
  #because ordinal variables are not yet supported we set these indicators to '0'
  #ordi <- numcats <- matrix(0,nrow=numG,ncol=P)
  cor.type <- lapply(1:numG,function(g){matrix(NA,ncol=P,nrow=P)})
  for(gg in 1:numG){
    for(p1 in 2:P){
      for(p2 in 1:(p1-1)){
        if(Ylevel[gg,p1]=="numeric" & Ylevel[gg,p2]=="numeric"){
          cor.type[[gg]][p1,p2]<-"product-moment"
        }else{
          if(Ylevel[gg,p1]=="dichotomous" & Ylevel[gg,p2]=="dichotomous"){
            cor.type[[gg]][p1,p2]<-"tetrachoric"
          }else{
            if(Ylevel[gg,p1]=="ordinal" & Ylevel[gg,p2]=="ordinal"){
              cor.type[[gg]][p1,p2]<-"polychoric"
            }else{
              if((Ylevel[gg,p1]=="ordinal" & Ylevel[gg,p2]=="numeric") |
                 (Ylevel[gg,p1]=="numeric" & Ylevel[gg,p2]=="ordinal")){
                cor.type[[gg]][p1,p2]<-"polyserial"
              }else{
                if((Ylevel[gg,p1]=="ordinal" & Ylevel[gg,p2]=="dichotomous")|
                   (Ylevel[gg,p1]=="dichotomous" & Ylevel[gg,p2]=="ordinal")){
                  cor.type[[gg]][p1,p2]<-"tetrachoric"
                }else{
                  if((Ylevel[gg,p1]=="numeric" & Ylevel[gg,p2]=="dichotomous")|
                     (Ylevel[gg,p1]=="dichotomous" & Ylevel[gg,p2]=="numeric")){
                    cor.type[[gg]][p1,p2]<-"biserial"
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  model_matrices <- lapply(seq_len(numG) , function(x) {
    model.matrix(formula, Y_groups[[x]])
  })

  correlate <- remove_predictors_helper(Y_groups = Y_groups, formula)

  YXlist <- lapply(1:length(model_matrices),function(g){
    list(as.matrix(correlate[[g]]),as.matrix(model_matrices[[g]]))
  })

  K <- ncol(YXlist[[1]][[2]])
  numcorr <- numG*P*(P-1)/2
  ngroups <- unlist(lapply(1:numG,function(g){nrow(YXlist[[g]][[1]])}))
  if(min(ngroups)<2*P+1){
    stop("sample size in each group must at least be equal to 2*P+1, with P the number of DVs.")
  }
  Ntot <- max(ngroups)
  Ygroups <- array(0,dim=c(numG,Ntot,P))
  Xgroups <- array(0,dim=c(numG,Ntot,K))
  XtXi <- array(0,dim=c(numG,K,K))
  BHat <- array(0,dim=c(numG,K,P))
  sdHat <- matrix(0,nrow=numG,ncol=P)
  CHat <- array(0,dim=c(numG,P,P))
  SumSq <- array(0,dim=c(numG,P,P))
  SumSqInv <- array(0,dim=c(numG,P,P))
  sdsd <- matrix(0,nrow=numG,ncol=P)

  for(g in 1:numG){
    Y_g <- YXlist[[g]][[1]]
    for(p in 1:P){
      if(ordi[g,p]==0){
        Y_g[,p] <- c(scale(Y_g[,p]))
      }
    }
    X_g <- YXlist[[g]][[2]]
    Ygroups[g,1:ngroups[g],] <- as.double(Y_g)
    #standardize data to get a more stable sampler for the correlations.
    tableX <- apply(X_g,2,table)
    catX <- unlist(lapply(1:length(tableX),function(xcol){
      length(tableX[[xcol]])
    }))
    if(sum(catX>1)){
      X_g[1:ngroups[g],which(catX>1)] <- apply(as.matrix(X_g[1:ngroups[g],which(catX>1)]),2,scale)
    }
    Xgroups[g,1:ngroups[g],] <- as.double(X_g)
    XtXi[g,,] <- as.double(solve(t(X_g)%*%X_g))
    BHat[g,,] <- as.double(XtXi[g,,]%*%t(X_g)%*%Y_g)
    SumSq[g,,] <- t(Y_g - X_g%*%BHat[g,,])%*%(Y_g - X_g%*%BHat[g,,])
    SumSqInv[g,,] <- solve(SumSq[g,,])
    Sigma_g <- SumSq[g,,]/ngroups[g]
    sdHat[g,] <- as.double(sqrt(diag(Sigma_g)))
    CHat[g,,] <- as.double(diag(1/sdHat[g,])%*%Sigma_g%*%diag(1/sdHat[g,]))
    #get rough estimate of posterior sd of the standard deviations (used for random walk sd)
    drawsSigma_g <- rWishart(1e2,df=ngroups[g],Sigma=SumSqInv[g,,])
    sdsd[g,] <- unlist(lapply(1:P,function(p){
      as.double(sd(sqrt(drawsSigma_g[p,p,])))
    }))
  }
  samsize0 <- iter
  gLiuSab <- array(as.double(0),dim=c(samsize0,numG,P))
  Njs <- matrix(as.double(ngroups),nrow=numG,ncol=1)

  # call Fortran subroutine for Gibbs sampling using noninformative improper priors
  # for regression coefficients, Jeffreys priors for standard deviations, and a proper
  # joint uniform prior for the correlation matrices.
  res <- .Fortran("estimate_bct_ordinal",
                  postZmean=matrix(as.double(0),nrow=numcorr,ncol=1),
                  postZcov=matrix(as.double(0),nrow=numcorr,ncol=numcorr),
                  P=as.integer(P),
                  numcorr=as.integer(numcorr),
                  K=as.integer(K),
                  numG=as.integer(numG),
                  BHat=BHat,
                  sdHat=sdsd,
                  CHat=CHat,
                  XtXi=XtXi,
                  samsize0=as.integer(samsize0),
                  burnin=as.integer(burnin),
                  Ntot=as.integer(Ntot),
                  Njs_in=Njs,
                  Xgroups=Xgroups,
                  Ygroups=Ygroups,
                  C_quantiles=array(as.double(0),dim=c(numG,P,P,3)),
                  sigma_quantiles=array(as.double(0),dim=c(numG,P,3)),
                  B_quantiles=array(as.double(0),dim=c(numG,K,P,3)),
                  BDrawsStore=array(as.double(0),dim=c(samsize0,numG,K,P)),
                  sigmaDrawsStore=array(as.double(0),dim=c(samsize0,numG,P)),
                  CDrawsStore=array(as.double(0),dim=c(samsize0,numG,P,P)),
                  sdMH=sdsd,
                  ordinal_in=ordi,
                  Cat_in=numcats,
                  maxCat=as.integer(max(numcats)),
                  gLiuSab=gLiuSab,
                  #seed=as.integer(sample.int(1e6,1)),
                  nuggetscale=as.double(nugget.scale),
                  priorchoice = as.integer(priorchoice)
                  #,WgroupsStore=array(as.double(0),dim=c(samsize0,numG,Ntot,P)),
                  #meanMatMeanStore = array(as.double(0),dim=c(samsize0,Ntot,P)),
                  #SigmaMatDrawStore = array(as.double(0),dim=c(samsize0,P,P)),
                  #CheckStore = array(as.double(1),dim=c(samsize0,numG,10,P,P))
  )

  varnames <- lapply(1:numG,function(g){
    names(correlate[[g]])
  })
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

  FmeansCovCorr <- lapply(1:numG,function(g){
    Fdraws_g <- FisherZ(t(matrix(unlist(lapply(1:samsize0,function(s){
      res$CDrawsStore[s,g,,][lower.tri(diag(P))]
    })),ncol=samsize0)))
    mean_g <- apply(Fdraws_g,2,mean)
    names(mean_g) <- corrnames[[g]][lower.tri(diag(P))]
    covm_g <- cov(Fdraws_g)
    ### DELETE THIS
    #covm_g <- diag(numcorr/numG)
    ### DELETE THIS
    return(list(mean_g,covm_g))
  })

  meansCovCorr <- lapply(1:numG,function(g){
    matcor_g <- unlist(lapply(1:(P-1),function(p2){
      unlist(lapply((p2+1):P,function(p1){
        mean(res$CDrawsStore[,g,p1,p2])
      }))
    }))
    names(matcor_g) <- corrnames[[g]][lower.tri(diag(P))]
    return(matcor_g)
  })

  meanN <- unlist(lapply(1:numG,function(g){
    FmeansCovCorr[[g]][[1]]
  }))
  covmN <- matrix(0,nrow=numcorr,ncol=numcorr)
  numcorrg <- numcorr/numG

  corrdraws <- lapply(1:numG,function(g){
    array_g <- res$CDrawsStore[,g,,]
    dimnames(array_g) <- list(NULL,varnames[[g]],varnames[[g]])
    return(array_g)
  })

  for(g in 1:numG){
    covmN[(g-1)*numcorrg+1:numcorrg,(g-1)*numcorrg+1:numcorrg] <- FmeansCovCorr[[g]][[2]]
    colnames(cor.type[[g]]) <- row.names(cor.type[[g]]) <- varnames[[g]]
  }

  # posterior estimates
  postestimates_correlations <- Reduce(rbind,
                                       lapply(1:numG,function(g){
                                         means <- meansCovCorr[[g]]
                                         medians <- res$C_quantiles[g,,,2][lower.tri(diag(P))]
                                         lb <- res$C_quantiles[g,,,1][lower.tri(diag(P))]
                                         ub <- res$C_quantiles[g,,,3][lower.tri(diag(P))]
                                         return(cbind(means,medians,lb,ub))
                                       }))
  colnames(postestimates_correlations) <- c("mean","median","2.5%","97.5%")

  cor_out <- list(meanF=meanN,covmF=covmN,correstimates=postestimates_correlations,
                  corrdraws=corrdraws,corrnames=corrnames,variables=varnames,
                  cor.type=cor.type,res=res,prior.cor=prior.cor)
  class(cor_out) <- "cor_test"

  return(cor_out)
}
