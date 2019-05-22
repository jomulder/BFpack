#' @importFrom pracma rref
#' @importFrom mvtnorm dmvnorm pmvnorm dmvt pmvt
#' @importFrom Matrix rankMatrix
#' @importFrom MASS ginv
#' @method BF glm
#' @export

source("./R/BFcorrelation.R")
source("./R/BFregression.R") #for make_Rrlist function

##TEST data
#multivariate lm data
ami_data <- read.table("http://static.lib.virginia.edu/statlab/materials/data/ami_data.DAT")
names(ami_data) <- c("TOT","AMI","GEN","AMT","PR","DIAP","QRS")
mlm1 <- lm(cbind(TOT, AMI) ~ GEN + AMT + PR + DIAP + QRS, data = ami_data)

constraints <- "TOT_with_AMI > AMI_with_AMI = TOT_with_TOT"

#glm data
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())

constraints <- "outcome2 = outcome3 > treatment2"


model <- glm.D93

var_names <- variable.names(model)
N <- length(model$residuals)
mean0 <- as.matrix(rep(0, length(var_names)))
covm0 <- vcov(model)
meanN <- coef(model)
covmN <- covm0*N

Gaussian_estimator <- function(mean0, covm0, meanN, covmN, constraints){

  # compute BFs and posterior probs using
  # prior mean en covmatrix mean0 and covm0
  # post mean en covmatrix meanN and covmN
  if(constraints=="exploratory"){ #IGNORING THIS PART FOR NOW> TALK TO JORIS TOMORROW
    # H0: corr = 0
    # H1: corr < 0
    # H2: corr < 0
    relfit <- matrix(c(dnorm(0,mean=meanN,sd=sqrt(diag(covmN))), #[Anton] Are these general or specific to correlations?
                       pnorm(0,mean=meanN,sd=sqrt(diag(covmN))),
                       1-pnorm(0,mean=meanN,sd=sqrt(diag(covmN)))),ncol=3)

    relcomp <- matrix(c(dnorm(0,mean=mean0,sd=sqrt(diag(covm0))),
                        pnorm(0,mean=mean0,sd=sqrt(diag(covm0))),
                        1-pnorm(0,mean=mean0,sd=sqrt(diag(covm0)))),ncol=3)
    BFtu <- relfit / relcomp
    PHP <- round(BFtu / apply(BFtu,1,sum),3)
    colnames(PHP) <- c("H0:corr=0","H1:corr<0","H2:corr>0")
    # get names of correlations
    corrmat <- diag(P)
    row.names(corrmat) <- colnames(corrmat) <- row.names(SumSquares[[1]])
    corr_names <- names(get_estimates(corrmat)$estimate)
    matrix_names <- matrix(corr_names,nrow=3)
    # equal correlations are at the opposite side of the vector
    corr_names <- matrix_names[lower.tri(matrix_names)]

    if(numgroup>1){
      matcorrpop <- matrix(0,nrow=length(corr_names),ncol=numgroup)
      for(c in 1:length(corr_names)){
        matcorrpop[c,] <- unlist(lapply(1:numgroup,function(pop){
          paste0(corr_names[c],"_gr",as.character(pop)) #or "_in_gr"
        }))
      }
      corr_names <- c(matcorrpop)
    }
    row.names(PHP) <- corr_names
    BFmatrix <- NULL

  }else{ #WORKING ON THIS PART
    # confirmatory tests based on input constraints

    #THIS IS ONLY ABOUT GETTING THE NAMES, SHOULD NOT BE INSIDE THIS PART OF THE FUNCTION
    # corrmat <- diag(P)
    # row.names(corrmat) <- colnames(corrmat) <- row.names(SumSquares[[1]])
    # corr_names <- names(get_estimates(corrmat)$estimate)
    # matrix_names <- matrix(corr_names,nrow=3)
    # # equal correlations are at the opposite side of the vector
    # corr_names <- c(matrix_names[lower.tri(matrix_names)],
    #                 t(matrix_names)[lower.tri(matrix_names)])
    # if(numgroup>1){
    #   matcorrpop <- matrix(0,nrow=length(corr_names),ncol=numgroup)
    #   for(c in 1:length(corr_names)){
    #     matcorrpop[c,] <- unlist(lapply(1:numgroup,function(pop){
    #       paste0(corr_names[c],"_gr",as.character(pop)) #or "_in_gr"
    #     }))
    #   }
    #   corr_names <- c(matcorrpop)
    # }

    #ADDITIONAL CORRELATION SPECIFIC (NAME EXTRACTION?)
    # select1 <- rep(1:numcorrgroup,numgroup) + rep((0:(numgroup-1))*2*numcorrgroup,each=numcorrgroup)
    # select2 <- rep(numcorrgroup+1:numcorrgroup,numgroup) + rep((0:(numgroup-1))*2*numcorrgroup,each=numcorrgroup)
    # #combine equivalent correlations, e.g., cor(Y1,Y2)=corr(Y2,Y1).
    # parse_hyp$hyp_mat <-
    #   cbind(parse_hyp$hyp_mat[,select1] + parse_hyp$hyp_mat[,select2],parse_hyp$hyp_mat[,numcorrgroup*2*numgroup+1])

    parse_hyp <- parse_hypothesis(var_names,constraints) #PARSE HYP REQUIRES COEF_NAMES

    #create coefficient with equality and order constraints
    RrList <- make_RrList2(parse_hyp) #[Anton]MAKE rRLIST 1 NOT WORKING AGAIN?
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]

    numhyp <- length(RrE)
    relcomp <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Gaussian_measures(mean0,covm0,RrE[[h]],RrO[[h]])
    })),nrow=2))
    relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
      Gaussian_measures(meanN,covmN,RrE[[h]],RrO[[h]])
    })),nrow=2))

    # get relative fit and complexity of complement hypothesis
    relcomp <- Gaussian_prob_Hc(mean0,covm0,relcomp,RrO) #[Anton]NOT WORKING!
    relfit <- Gaussian_prob_Hc(meanN,covmN,relfit,RrO)

    Hnames <- c(unlist(lapply(1:numhyp,function(h){paste0("H",as.character(h))})),"Hc")
    row.names(relcomp) <- Hnames
    row.names(relfit) <- Hnames


    #[Anton]QUESTION HERE! Should this be inside or outside the function? I.e., could all functions use the same
    #summarizing process (approximately below) after computing relcomp and relfit?

    # the BF for the complement hypothesis vs Hu needs to be computed.
    BFtu <- c(apply(relfit / relcomp, 1, prod))
    # Check input of prior probabilies
    if(!(priorprob == "default" || (length(priorprob)==nrow(relfit) && min(priorprob)>0) )){
      stop("'probprob' must be a vector of positive values or set to 'default'.")
    }
    # Change prior probs in case of default setting
    if(priorprob=="default"){priorprobs <- rep(1,length(BFtu))}
    PHP <- round(BFtu*priorprobs / sum(BFtu*priorprobs),3)
    BFmatrix <- matrix(rep(BFtu,length(BFtu)),ncol=length(BFtu))/
      t(matrix(rep(BFtu,length(BFtu)),ncol=length(BFtu)))
    row.names(BFmatrix) <- Hnames
    colnames(BFmatrix) <- Hnames
  }

  return(list(BFtu=BFtu,PHP=PHP,BFmatrix=BFmatrix,relfit=relfit,relcomp=relcomp,
              SumSquares=SumSquares,n=n,constraints=constraints,nu0=nu0,mean0=mean0,
              covm0=covm0,priorprob=priorprob,corr_names=corr_names))
}
