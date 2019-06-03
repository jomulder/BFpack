###Default method for BFpack


#' @method BF default
#' @export
BF.default <- function(x,
                       hypothesis = NULL,
                       prior = NULL,
                       ...){

  #Extract model data
  names_coef <- variable.names(x)
  N <- nobs(x)
  mean0 <- as.matrix(rep(0, length(names_coef)))
  covm0 <- vcov(x) * N
  meanN <- coef(x)
  covmN <- vcov(x)

  # compute BFs and posterior probs using
  # prior mean and covmatrix mean0 and covm0
  # post mean and covmatrix meanN and covmN
    # H0: corr = 0
    # H1: corr < 0
    # H2: corr < 0
    relfit <- matrix(c(dnorm(0,mean=meanN,sd=sqrt(diag(covmN))), #[Anton] Are these relfit/relcomp computations general or specific to correlations? [Joris] This is general. So it also works for these parameters.
                       pnorm(0,mean=meanN,sd=sqrt(diag(covmN))),
                       1-pnorm(0,mean=meanN,sd=sqrt(diag(covmN)))),ncol=3)

    relcomp <- matrix(c(dnorm(0,mean=mean0,sd=sqrt(diag(covm0))),
                        pnorm(0,mean=mean0,sd=sqrt(diag(covm0))),
                        1-pnorm(0,mean=mean0,sd=sqrt(diag(covm0)))),ncol=3)
    BFtu_exploratory <- relfit / relcomp
    PHP_exploratory <- round(BFtu_exploratory / apply(BFtu_exploratory,1,sum),3)
    colnames(PHP_exploratory) <- c("p(=0)","Pr(<0)","Pr(>0)")
    row.names(PHP_exploratory) <- names_coef

  if(!is.null(hypothesis)){
    # confirmatory tests based on input constraints

    parse_hyp <- parse_hypothesis(names_coef, hypothesis)

    #create coefficient with equality and order constraints
    RrList <- make_RrList2(parse_hyp)
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]

    #get relative fit and complexity of hypotheses
    numhyp <- length(RrE)
    relcomp <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Gaussian_measures(mean0,covm0,RrE[[h]],RrO[[h]])
    })),nrow=2))

    relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Gaussian_measures(meanN,covmN,RrE[[h]],RrO[[h]])
    })),nrow=2))

    # get relative fit and complexity of complement hypothesis
    relcomp <- Gaussian_prob_Hc(mean0,covm0,relcomp, RrO, RrO) #Note that input is a bit strange here, Gaussian_prob_Hc needs fixing
    relfit <- Gaussian_prob_Hc(meanN,covmN,relfit,RrO, RrO)

    Hnames <- c(unlist(lapply(1:numhyp,function(h){paste0("H",as.character(h))})),"Hc")
    row.names(relcomp) <- row.names(relfit) <- Hnames
    colnames(relcomp) <- c("c_E", "c_0")
    colnames(relfit) <- c("f_E", "f_0")


    # the BF for the complement hypothesis vs Hu needs to be computed.
    BFtu_confirmatory <- c(apply(relfit / relcomp, 1, prod))
    # Check input of prior probabilies
    if(!(is.null(prior) || (length(prior)==nrow(relfit) && min(prior)>0) )){
      stop("'prior' must be a vector of positive values or set to 'NULL'.")
    }
    # Change prior probs in case of default setting
    if(is.null(prior)){prior <- rep(1,length(BFtu_confirmatory))}

    PHP_confirmatory <- round(BFtu_confirmatory*prior / sum(BFtu_confirmatory*prior),3)
    BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory))/
      t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
    row.names(BFmatrix_confirmatory) <- Hnames
    colnames(BFmatrix_confirmatory) <- Hnames
  }

    out <- list(
      BFtu_exploratory=BFtu_exploratory,
      PHP_exploratory=PHP_exploratory,
      BFtu_confirmatory=BFtu_confirmatory,
      PHP_confirmatory=PHP_confirmatory,
      BFmatrix_confirmatory=BFmatrix_confirmatory,
      relative_fit=relfit,
      relative_complexity=relcomp,
      model=x,
      estimates=meanN,
      constraints=hypothesis,
      priorprob=prior)

    class(out) <- "BF"

    out

}

