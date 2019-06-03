#'
#'
#'
#'
#' #' @importFrom pracma rref
#' #' @importFrom mvtnorm dmvnorm pmvnorm rmvnorm dmvt pmvt
#' #' @importFrom Matrix rankMatrix
#' #' @importFrom MCMCpack rinvgamma
#' #' @importFrom MASS ginv
#' #' @method BF mlm
#' #' @export
#'
#' BF.mlm <- function(x,
#'                   hypothesis = NULL,
#'                   prior = NULL,
#'                   covariates = NULL,
#'                   parameter = NULL,
#'                   ...){
#'
#'   if(is.null(hypothesis)){
#'     constraints <- "exploratory"
#'   } else {
#'     constraints <- hypothesis
#'   }
#'   if(is.null(prior)){
#'     priorprob <- "default"
#'   } else {
#'     priorprob <- prior
#'   }
#'   if(is.null(parameter)){
#'     parametertest <- "regression"
#'   }else if(parameter=="correlation"){
#'     parametertest <- parameter
#'   }else{parametertest <- "regression"}
#'
#'
#'   # default BF on location parameters in a univarite normal linear model
#'   # Note that it is recommended that the fitten model is based on standardized covariates.
#'   if(!is.matrix(x$coefficients)){
#'     P <- 1
#'     N <- length(x$residuals)
#'     K <- length(x$coefficients) # dimension of predictors
#'     dummyX <- rep(F,K)
#'     names(dummyX) <- names(x$coefficients)
#'   } else {
#'     P <- ncol(x$residuals)
#'     N <- nrow(x$residuals)
#'     K <- length(x$coefficients)/P # dimension of predictors per dependent variable
#'     dummyX <- rep(F,K)
#'     names(dummyX) <- row.names(x$coefficients)
#'   }
#'
#'   if(is.null(covariates)){
#'     noncovs <- 1:K
#'   } else { # covariates must be a vector of integers denoting which predictor variables
#'     # are not grouping variables.
#'     noncovs <- (1:K)[-covariates]
#'   }
#'
#'   Xmat <- model.matrix(x)
#'   Ymat <- model.matrix(x)%*%x$coefficients + x$residuals
#'
#'   if(parametertest=="regression"){
#'     #Check which are dummy variables corresponding to (adjusted) mean parameters
#'     for(k in noncovs){
#'       uniquek <- sort(unique(Xmat[,k]))
#'       # if(length(uniquek)==2){
#'       #   if(uniquek[1]==0 && uniquek[2]==1){dummyX[k]<-T} #(adjusted) mean
#'       # }else{
#'       #   if(length(uniquek)==1){dummyX[k]<-T} #intercept parameter
#'       # }
#'       if(length(uniquek)<=2){dummyX[k]<-T} #group index of intercept
#'     }
#'     #number of groups on variations of dummy combinations
#'     groupcode <- as.matrix(unique(Xmat[,dummyX]))
#'     rownames(groupcode) <- unlist(lapply(1:nrow(groupcode),function(r){
#'       paste0("groupcode",r)
#'     }))
#'
#'     J <- nrow(groupcode)
#'     if(J==0){ #then no intercept
#'       J <- 1
#'       nointercept <- T
#'     }
#'     # group membership of each observation
#'     dvec <- unlist(lapply(1:N,function(i){
#'       which(rowSums(abs(t(matrix(rep(Xmat[i,dummyX],J),ncol=J)) - groupcode))==0)
#'     }))
#'     Nj <- c(table(dvec))
#'     #set minimal fractions for each group
#'     bj <- ((P+K)/J)/Nj
#'
#'     #Compute sufficient statistics for all groups
#'     tXXj <- lapply(1:J,function(j){
#'       if(Nj[j]==1){
#'         Xmat[dvec==j,]%*%t(Xmat[dvec==j,])
#'       }else t(Xmat[dvec==j,])%*%Xmat[dvec==j,]
#'     })
#'     tXXj_b <- lapply(1:J,function(j){
#'       tXXj[[j]]*bj[j]
#'     })
#'     tXYj <- lapply(1:J,function(j){
#'       if(Nj[j]==1){
#'         as.matrix(Xmat[dvec==j,]*Ymat[dvec==j,])
#'       } else {t(Xmat[dvec==j,])%*%Ymat[dvec==j,]}
#'     })
#'     tXYj_b <- lapply(1:J,function(j){
#'       tXYj[[j]]*bj[j]
#'     })
#'     tYYj <- lapply(1:J,function(j){
#'       t(Ymat[dvec==j,])%*%Ymat[dvec==j,]
#'     })
#'     tYYj_b <- lapply(1:J,function(j){
#'       tYYj[[j]]*bj[j]
#'     })
#'     tXX <- Reduce("+",tXXj)
#'     tXXi <- solve(tXX)
#'     tXY <- Reduce("+",tXYj)
#'     tYY <- Reduce("+",tYYj)
#'     tXX_b <- Reduce("+",tXXj_b)
#'     tXXi_b <- solve(tXX_b)
#'     tXY_b <- Reduce("+",tXYj_b)
#'     tYY_b <- Reduce("+",tYYj_b)
#'     BetaHat <- solve(tXX)%*%tXY          # same as x$coefficients
#'     S <- tYY - t(tXY)%*%solve(tXX)%*%tXY # same as sum((x$residuals)**2)
#'     # sufficient statistics based on fraction of the data
#'     BetaHat_b <- solve(tXX_b)%*%tXY_b
#'     S_b <- tYY_b - t(tXY_b)%*%solve(tXX_b)%*%tXY_b
#'
#'
#'     # BF computation for exploratory analysis of separate parameters
#'     if(P==1){
#'       names_coef <- names(x$coefficients)
#'     }else{
#'       names_coef1 <- names(x$coefficients[,1])
#'       names_coef2 <- names(x$coefficients[1,])
#'       names_coef <- unlist(lapply(1:P,function(p){
#'         lapply(1:K,function(k){
#'           paste(names_coef1[k],".",names_coef2[p],sep="")
#'         })
#'       }))
#'     }
#'
#'     # prior hyperparameters
#'     df0 <- 1 # should be the same as sum(rep(bj,times=Nj))-K-P+1
#'     Scale0 <- kronecker(S_b,tXXi_b)
#'     mean0 <- as.matrix(rep(0,K*P))
#'     # posterior hyperparameters
#'     dfN <- N-K-P+1
#'     ScaleN <- kronecker(S,tXXi)/(N-K-P+1) # off-diagonal elements have no meaning
#'     meanN <- as.matrix(c(BetaHat))
#'
#'     # Hypotheses for exploratory test
#'     # H0: beta = 0
#'     # H1: beta < 0
#'     # H2: beta < 0
#'     relfit <- t(matrix(unlist(lapply(1:(K*P),function(k){
#'       c(dt((0-meanN[k,1])/sqrt(ScaleN[k,k]),df=dfN)/sqrt(ScaleN[k,k]),
#'         pt((0-meanN[k,1])/sqrt(ScaleN[k,k]),df=dfN,lower.tail = TRUE),
#'         pt((0-meanN[k,1])/sqrt(ScaleN[k,k]),df=dfN,lower.tail = FALSE))
#'     })),nrow=3))
#'     relcomp <- t(matrix(unlist(lapply(1:(K*P),function(k){
#'       c(dt((0-mean0[k,1])/sqrt(Scale0[k,k]),df=df0)/sqrt(Scale0[k,k]),
#'         pt((0-mean0[k,1])/sqrt(Scale0[k,k]),df=df0,lower.tail = TRUE),
#'         pt((0-mean0[k,1])/sqrt(Scale0[k,k]),df=df0,lower.tail = FALSE))
#'     })),nrow=3))
#'     colnames(relfit) <- colnames(relcomp) <- c("p(=0)","Pr(<0)","Pr(>0)")
#'     row.names(relcomp) <- row.names(relfit) <- names_coef
#'
#'     BFtu_exploratory <- relfit / relcomp
#'     colnames(BFtu_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
#'     PHP_exploratory <- BFtu_exploratory / apply(BFtu_exploratory,1,sum)
#'
#'     if(constraints!="exploratory"){
#'       #read constraints
#'       if(P==1){
#'         names_coef <- names(x$coefficients)
#'       }else{
#'         names_coef1 <- names(x$coefficients[,1])
#'         names_coef2 <- names(x$coefficients[1,])
#'         names_coef <- unlist(lapply(1:P,function(p){
#'           lapply(1:K,function(k){
#'             paste(names_coef1[k],".",names_coef2[p],sep="")
#'           })
#'         }))
#'       }
#'       # translate named constraints to matrices with coefficients for constraints
#'       parse_hyp <- parse_hypothesis(names_coef,constraints)
#'       RrList <- make_RrList2(parse_hyp)
#'       RrE <- RrList[[1]]
#'       RrO <- RrList[[2]]
#'
#'       RrStack <- rbind(do.call(rbind,RrE),do.call(rbind,RrO))
#'       RStack <- RrStack[,-(K+1)]
#'       rStack <- RrStack[,(K+1)]
#'
#'       # check if a common boundary exists for prior location under all constrained hypotheses
#'       if(nrow(RrStack) > 1){
#'         rref_ei <- rref(RrStack)
#'         nonzero <- RrStack[,K+1]!=0
#'         if(max(nonzero)>0){
#'           row1 <- max(which(nonzero==T))
#'           if(sum(abs(RrStack[row1,1:K]))==0){
#'             stop("No common boundary point for prior location. Conflicting constraints.")
#'           }
#'         }
#'       }
#'
#'       if(P==1){ # then the unconstrained prior and posterior have multivariate Student t distributions
#'
#'         # prior hyperparameters
#'         df0 <- 1 # should be the same as sum(rep(bj,times=Nj))-K-P+1
#'         Scale0 <- kronecker(S_b,tXXi_b)
#'         mean0 <- as.matrix(rep(0,K*P))
#'         # posterior hyperparameters
#'         dfN <- N-K-P+1
#'         ScaleN <- kronecker(S,tXXi)/(N-K-P+1) # off-diagonal elements have no meaning
#'         meanN <- as.matrix(c(BetaHat))
#'
#'         #number of hypotheses that are specified
#'         numhyp <- length(RrO)
#'
#'         relcomp <- t(matrix(unlist(lapply(1:numhyp,function(h){
#'           Student_measures(mean0,Scale0,df0,RrE[[h]],RrO[[h]])
#'         })),nrow=2))
#'
#'         relfit <- t(matrix(unlist(lapply(1:numhyp,function(h){
#'           Student_measures(meanN,ScaleN,dfN,RrE[[h]],RrO[[h]])
#'         })),nrow=2))
#'
#'         # Compute relative fit/complexity for the complement hypothesis
#'         relfit <- Student_prob_Hc(meanN,scaleN,dfN,relfit,constraints)
#'         relcomp <- Student_prob_Hc(mean0,scale0,df0,relcomp,constraints)
#'         row.names(relcomp)[1:numhyp] <- parse_hyp$original_hypothesis
#'         row.names(relfit)[1:numhyp] <- parse_hyp$original_hypothesis
#'         colnames(relcomp) <- c("c_E","c_O")
#'         colnames(relfit) <- c("f_E","f_O")
#'
#'       }else{
#'
#'         #number of hypotheses that are specified
#'         numhyp <- length(RrO)
#'         Mean0 <- matrix(0,nrow=K,ncol=P)
#'
#'         relmeasunlist <- unlist(lapply(1:numhyp,function(h){
#'           # Check whether the constraints are on a single row or column, if so
#'           # use the analytic expression, else using a Monte Carlo estimate.
#'           RrStack <- rbind(RrE[[h]],RrO[[h]])
#'           Rcheck <- Reduce("+",lapply(1:nrow(RrStack),function(row1){
#'             abs(matrix(RrStack[row1,-(K*P+1)],nrow=K))
#'           }))
#'           RcheckRow <- apply(Rcheck,1,sum)
#'           RcheckCol <- apply(Rcheck,2,sum)
#'           if(sum(RcheckRow!=0)==1){ # use multivariate Student distributions
#'             K1 <- which(RcheckRow!=0)
#'             # posterior hyperparameters
#'             dfN <- N-K-P+1
#'             ScaleN <- S*tXXi[K1,K1]/(N-K-P+1) # off-diagonal elements have no meaning
#'             meanN <- as.matrix(c(BetaHat[K1,]))
#'             # exclude inactive rows
#'             if(is.null(RrE[[h]])){RrE_h=NULL
#'             }else{
#'               if(nrow(RrE[[h]])==1){
#'                 RrE_h <- t(as.matrix(RrE[[h]][,c((0:(P-1))*K+K1,P*K+1)]))
#'               }else{
#'                 RrE_h <- RrE[[h]][,c((0:(P-1))*K+K1,P*K+1)]
#'               }
#'             }
#'             if(is.null(RrO[[h]])){RrO_h=NULL
#'             }else{
#'               if(nrow(RrO[[h]])==1){
#'                 RrO_h <- t(as.matrix(RrO[[h]][,c((0:(P-1))*K+K1,P*K+1)]))
#'               }else{
#'                 RrO_h <- RrO[[h]][,c((0:(P-1))*K+K1,P*K+1)]
#'               }
#'             }
#'             # prior hyperparameters
#'             df0 <- 1 # should be the same as sum(rep(bj,times=Nj))-K-P+1
#'             Scale0 <- S_b*tXXi_b[K1,K1]
#'             mean0 <- Mean0[K1,]
#'             # compute relative measures of fit and complexity
#'             relcomp_h <- Student_measures(mean0,Scale0,df0,RrE_h,RrO_h)
#'             relfit_h <- Student_measures(meanN,ScaleN,dfN,RrE_h,RrO_h)
#'
#'           }else if(sum(RcheckCol!=0)==1){ # use multivariate Student distributions
#'             P1 <- which(RcheckCol!=0)
#'             # prior hyperparameters
#'             df0 <- 1 # should be the same as sum(rep(bj,times=Nj))-K-P+1
#'             Scale0 <- S_b[P1,P1]*tXXi_b
#'             mean0 <- Mean0[,P1]
#'             # posterior hyperparameters
#'             dfN <- N-K-P+1
#'             ScaleN <- S[P1,P1]*tXXi/(N-K-P+1) # off-diagonal elements have no meaning
#'             meanN <- as.matrix(c(BetaHat[,P1]))
#'             # exclude inactive rows
#'             if(is.null(RrE[[h]])){RrE_h=NULL
#'             }else{
#'               if(nrow(RrE[[h]])==1){
#'                 RrE_h <- t(as.matrix(RrE[[h]][,c((P1-1)*K+1:K,P*K+1)]))
#'               }else{
#'                 RrE_h <- RrE[[h]][,c((P1-1)*K+1:K,P*K+1)]
#'               }
#'             }
#'             if(is.null(RrO[[h]])){RrO_h=NULL
#'             }else{
#'               if(nrow(RrO[[h]])==1){
#'                 RrO_h <- t(as.matrix(RrO[[h]][,c((P1-1)*K+1:K,P*K+1)]))
#'               }else{
#'                 RrO_h <- RrO[[h]][,c((P1-1)*K+1:K,P*K+1)]
#'               }
#'             }
#'             # compute relative measures of fit and complexity
#'             relcomp_h <- Student_measures(mean0,Scale0,df0,RrE_h,RrO_h)
#'             relfit_h <- Student_measures(meanN,ScaleN,dfN,RrE_h,RrO_h)
#'
#'           }else{ #use Matrix-Student distributions with Monte Carlo estimate
#'
#'             df0 <- 1
#'             dfN <- N-K-P+1
#'             relfit_h <- MatrixStudent_measures(BetaHat,S,tXXi,dfN,RrE[[h]],RrO[[h]],MCdraws=1e4)
#'             relcomp_h <- MatrixStudent_measures(Mean0,S_b,tXXi_b,df0,RrE[[h]],RrO[[h]],MCdraws=1e4)
#'           }
#'           return(list(relfit_h,relcomp_h))
#'         }))
#'
#'         relfit <- t(matrix(unlist(relmeasunlist)[rep((0:(numhyp-1))*4,each=2)+rep(1:2,numhyp)],nrow=2))
#'         row.names(relfit) <- parse_hyp$original_hypothesis
#'         colnames(relfit) <- c("f_E","f_O")
#'         relcomp <- t(matrix(unlist(relmeasunlist)[rep((0:(numhyp-1))*4,each=2)+rep(3:4,numhyp)],nrow=2))
#'         row.names(relcomp) <- parse_hyp$original_hypothesis
#'         colnames(relcomp) <- c("c_E","c_O")
#'
#'         # Compute relative fit/complexity for the complement hypothesis
#'         relfit <- MatrixStudent_prob_Hc(BetaHat,S,tXXi,N-K-P+1,as.matrix(relfit),RrO)
#'         relcomp <- MatrixStudent_prob_Hc(Mean0,S_b,tXXi_b,1,as.matrix(relcomp),RrO)
#'       }
#'
#'       # the BF for the complement hypothesis vs Hu needs to be computed.
#'       BFtu_confirmatory <- c(apply(relfit / relcomp, 1, prod))
#'       # Check input of prior probabilies
#'       if(!(priorprob == "default" || (length(priorprob)==nrow(relfit) && min(priorprob)>0) )){
#'         stop("'probprob' must be a vector of positive values or set to 'default'.")
#'       }
#'       # Change prior probs in case of default setting
#'       if(priorprob=="default"){
#'         priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
#'       }else{
#'         priorprobs <- priorprobs/sum(priorprobs)
#'       }
#'       names(priorprobs) <- names(BFtu_confirmatory)
#'       PHP_confirmatory <- BFtu_confirmatory*priorprobs / sum(BFtu_confirmatory*priorprobs)
#'       # names(PHP_confirmatory) <- unlist(lapply(1:length(parse_hyp$original_hypothesis),function(hyp){
#'       #   paste0("Pr(",parse_hyp$original_hypothesis,")")
#'       # }))
#'       BFmatrix_confirmatory <- matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory))/
#'         t(matrix(rep(BFtu_confirmatory,length(BFtu_confirmatory)),ncol=length(BFtu_confirmatory)))
#'       row.names(BFmatrix_confirmatory) <- colnames(BFmatrix_confirmatory) <- names(BFtu_confirmatory)
#'     }else{
#'       BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relfit <-
#'         relcomp <- NULL
#'     }
#'   }else{ #perform correlation tests
#'
#'     # corrmat <- diag(P)
#'     # row.names(corrmat) <- colnames(corrmat) <- colnames(x$residuals)
#'     # corr_names <- names(get_estimates(corrmat)$estimate)
#'     # matrix_names <- matrix(corr_names,nrow=3)
#'     # # equal correlations are at the opposite side of the vector
#'     # corr_names <- c(matrix_names[lower.tri(matrix_names)],
#'     #                 t(matrix_names)[lower.tri(matrix_names)])
#'     # params_in_hyp1 <- params_in_hyp(constraints)
#'     # if(abs(unlist(lapply(1:length(params_in_hyp1),function(par){
#'     #   sum(grepl(params_in_hyp1[par],corr_names))
#'     # }))-rep(1,length(params_in_hyp1)))==0){
#'     #   ngroups <- 1
#'     # }else{ #check if labels of correlations in hypotheses match with
#'     #   #combined labels of correlations and dummy group variables
#'     #
#'     #   #Check which are dummy variables corresponding to (adjusted) mean parameters
#'     #   for(k in 1:K){
#'     #     uniquek <- sort(unique(Xmat[,k]))
#'     #     if(length(uniquek)<=2){dummyX[k]<-T} #group index of intercept
#'     #   }
#'     #   if(sum(dummyX)==0){
#'     #     # there appears to be no dummy group variables
#'     #     stop("Hypotheses appear to be formulated on correlations with incorrect labels. Labels for correlations should be of the form 'y1_with_y2' or 'y1_with_y2_in_x1'.")
#'     #   }
#'     #   # create correlation labels
#'     #   corr_names_all <- unlist(lapply(which(dummyX==T),function(k){
#'     #     if(dummyX[k]){
#'     #       unlist(lapply(1:length(corr_names),function(naam){
#'     #         paste0(corr_names[naam],"_in_",names(dummyX[k]))
#'     #       }))
#'     #     }
#'     #   }))
#'     #   names(corr_names_all) <- NULL
#'     #   unlist(lapply(1:length(params_in_hyp1),function(par){
#'     #     sum(grepl(params_in_hyp1[par],corr_names_all))
#'     #   }))
#'     #
#'     #   ## ABOVE SHOULD CHECK IF EXACT NAMES MATCH NOT PART OF NAMES!!!!!
#'     #
#'     #
#'     #
#'     #
#'     #
#'     # }
#'
#'
#'
#'
#'
#'
#'
#'     }
#'
#'     if(constraints!="exploratory"){
#'       # check if there are group specific correlations
#'
#'     }
#'     if(numgroup>1){
#'       matcorrpop <- matrix(0,nrow=length(corr_names),ncol=numgroup)
#'       for(c in 1:length(corr_names)){
#'         matcorrpop[c,] <- unlist(lapply(1:numgroup,function(pop){
#'           paste0(corr_names[c],"_gr",as.character(pop)) #or "_in_gr"
#'         }))
#'       }
#'       corr_names <- c(matcorrpop)
#'     }
#'
#'     parse_hyp <- parse_hypothesis(corr_names,constraints)
#'     select1 <- rep(1:numcorrgroup,numgroup) + rep((0:(numgroup-1))*2*numcorrgroup,each=numcorrgroup)
#'     select2 <- rep(numcorrgroup+1:numcorrgroup,numgroup) + rep((0:(numgroup-1))*2*numcorrgroup,each=numcorrgroup)
#'     #combine equivalent correlations, e.g., cor(Y1,Y2)=corr(Y2,Y1).
#'     parse_hyp$hyp_mat <-
#'       cbind(parse_hyp$hyp_mat[,select1] + parse_hyp$hyp_mat[,select2],parse_hyp$hyp_mat[,numcorrgroup*2*numgroup+1])
#'     #create coefficient with equality and order constraints
#'     RrList <- make_RrList(parse_hyp)
#'     RrE <- RrList[[1]]
#'     RrO <- RrList[[2]]
#'
#'
#'
#'
#'   }
#'
#'   BFlm_out <- list(
#'     BFtu_exploratory=BFtu_exploratory,
#'     PHP_exploratory=PHP_exploratory,
#'     BFtu_confirmatory=BFtu_confirmatory,
#'     PHP_confirmatory=PHP_confirmatory,
#'     BFmatrix_confirmatory=BFmatrix_confirmatory,
#'     relative_fit=relfit,
#'     relative_complexity=relcomp,
#'     model=x,
#'     estimates=x$coefficients,
#'     P=P,
#'     ngroups=Nj,
#'     tXXj=tXXj,
#'     tXYj=tXYj,
#'     tYYj=tYYj,
#'     constraints=constraints,
#'     priorprob=priorprob,
#'     groupcode=groupcode,
#'     dummyX=dummyX)
#'
#'   class(BFlm_out) <- "BF"
#'
#'   return(BFlm_out)
#' }
#'
#'
#'
#' params_in_hyp <- function(hyp){
#'   params_in_hyp <- trimws(unique(strsplit(hyp, split = "[ =<>,\\(\\);&\\*+-]+", perl = TRUE)[[1]]))
#'   params_in_hyp <- params_in_hyp[!sapply(params_in_hyp, grepl, pattern = "^[0-9]*\\.?[0-9]+$")]
#'   params_in_hyp[grepl("^[a-zA-Z]", params_in_hyp)]
#' }
#'
#'
