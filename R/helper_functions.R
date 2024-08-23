check_vcov <- function(x){
  if (!isTRUE(all.equal(x, t(x))) || any(diag(x) < 0)){
    saveRDS(x, "c:/git_repositories/BFpack/erordump.RData")
    stop(sQuote("sigma"), " is not a covariance matrix")
  }
}

process.prior.hyp.explo <- function(prior_hyp_explo, model){

  if(sum(class(model)=="bartlett_htest")>0){
    if(is.null(prior_hyp_explo)){
      prior_hyp_explo <- list(rep(.5,2))
    }else{
      if(!is.list(prior_hyp_explo)){
        prior_hyp_explo <- list(prior_hyp_explo)
      }else{
        prior_hyp_explo <- list(prior_hyp_explo[[1]])
      }
      if(length(prior_hyp_explo[[1]])!=2){
        stop("For an object of class 'bartlett_htest', the argument 'prior_hyp_explo' should be a vector
             of length 2 of the prior probabilities for the hypotheses of homogeneity of variances and an
             unconstrained alternative. Or use the default ('NULL') which implies equal prior probabilities.")
      }
    }
  }else{
    if(is.null(prior_hyp_explo)){ # then equal prior probabilities
      if(sum(class(model)=="aov")>0){
        prior_hyp_explo <- list(c(.5,.25,.25),rep(1,2),rep(1,2))
      }else{
        prior_hyp_explo <- list(c(.5,.25,.25))
      }
    }else{
      if(!is.list(prior_hyp_explo)){#then no error if it is a vector of length 3 for testing individual parameters
        if(length(prior_hyp_explo)==3){
          if(sum(class(model)=="aov")>0){
            prior_hyp_explo <- list(prior_hyp_explo,rep(1,2),rep(1,2))
          }else{
            prior_hyp_explo <- list(prior_hyp_explo)
          }
        }else{
          stop("The argument 'prior_hyp_explo' must be a vector of length three specifying the
             prior probabilities of a zero, negative, or positive parameter. In the case of an object of
             class 'aov' then 'prior_hyp_explo' can be a list of three elements containing the three prior
             probabilities of the exploratory test of the parameters, the two prior
             probabilities for testing the main effects, and the two prior probabilities for testing the
             interaction effects (check the documentation: ?BF).")
        }
      }else{
        if(length(prior_hyp_explo)==1){
          if(is.null(prior_hyp_explo[[1]])){
            prior_hyp_explo <- list(c(.5,.25,.25),rep(1,2),rep(1,2))
          }else{
            if(length(prior_hyp_explo[[1]]!=3)){
              stop("Specify three prior probabilities for the exploratory test of a zero, negative, or positive
                 effect or use the default 'NULL' implying equal prior probabilities (check the documentation: ?BF).")
            }
            prior_hyp_explo <- lapply(prior_hyp_explo,function(x){x/sum(x)})
          }
        }else{
          if(sum(class(model)=="aov")>0){

            if(length(prior_hyp_explo)==2){
              if(is.null(prior_hyp_explo[[1]])){
                prior_hyp_explo[[1]] <- c(.5,.25,.25)
              }else{
                if(length(prior_hyp_explo[[1]]!=3)){
                  stop("Specify three prior probabilities for the exploratory test of a zero, negative, or positive
                 effect or use the default 'NULL' implying equal prior probabilities (check the documentation: ?BF).")
                }
              }
              if(is.null(prior_hyp_explo[[2]])){
                prior_hyp_explo[[2]] <- rep(1,2)
              }else{
                if(length(prior_hyp_explo[[2]]!=2)){
                  stop("Specify two prior probabilities for the exploratory testing of main effects
                   or use the default 'NULL' implying equal prior probabilities
                   (check the documentation: ?BF).")
                }
              }
              prior_hyp_explo[[3]] <- rep(1,2)
            }
            if(length(prior_hyp_explo)==3){
              if(is.null(prior_hyp_explo[[1]])){
                prior_hyp_explo[[1]] <- c(.5,.25,.25)
              }else{
                if(length(prior_hyp_explo[[1]])!=3){
                  stop("Specify three prior probabilities for the exploratory test of a zero, negative, or positive
                 effect or use the default 'NULL' implying equal prior probabilities
                   (check the documentation: ?BF).")
                }
              }
              if(is.null(prior_hyp_explo[[2]])){
                prior_hyp_explo[[2]] <- rep(1,2)
              }else{
                if(length(prior_hyp_explo[[2]])!=2){
                  stop("Specify two prior probabilities for the exploratory testing of main effects
                   or use the default 'NULL' implying equal prior probabilities
                   (check the documentation: ?BF).")
                }
              }
              if(is.null(prior_hyp_explo[[3]])){
                prior_hyp_explo[[3]] <- rep(1,2)
              }else{
                if(length(prior_hyp_explo[[3]])!=2){
                  stop("Specify two prior probabilities for the exploratory testing of interaction effects
                   or use the default 'NULL' implying equal prior probabilities
                   (check the documentation: ?BF).")
                }
              }
            }
            if(length(prior_hyp_explo)>3){
              stop("Use a list of three vectors or use the default 'NULL' implying equal prior probabilities
               (check the documentation: ?BF).")
            }
          }else{
            if(is.null(prior_hyp_explo[[1]])){
              prior_hyp_explo <- list(c(.5,.25,.25))
            }else{
              if(length(prior_hyp_explo[[1]])!=3){
                stop("Specify three prior probabilities for the exploratory test of a zero, negative, or positive
                 effect or use the default 'NULL' implying equal prior probabilities (check the documentation: ?BF).")
              }
              prior_hyp_explo <- list(prior_hyp_explo[[1]]/sum(prior_hyp_explo[[1]]))
            }
          }
        }
      }
    }
    prior_hyp_explo <- lapply(prior_hyp_explo,function(x)x/sum(x))
  }

  return(prior_hyp_explo)

}

# from the output of the constraints in 'parse_hypothesis' create lists for the equality and order matrices
make_RrList <- function(parse_hyp){
  numhyp <- length(parse_hyp$hyp_mat)
  RrE <- lapply(1:numhyp,function(h){
    qE <- parse_hyp$n_constraints[h*2-1]
    if(qE==1){
      RrE_h <- t(as.matrix(parse_hyp$hyp_mat[[h]][1:qE,]))
    }else if(qE>1){
      RrE_h <- parse_hyp$hyp_mat[[h]][1:qE,]
    }else {RrE_h=NULL}
    RrE_h
  })
  RrO <- lapply(1:numhyp,function(h){
    qE <- parse_hyp$n_constraints[h*2-1]
    qO <- parse_hyp$n_constraints[h*2]
    if(qO==1){
      RrO_h <- t(as.matrix(parse_hyp$hyp_mat[[h]][qE+1:qO,]))
    }else if(qO>1){
      RrO_h <- parse_hyp$hyp_mat[[h]][qE+1:qO,]
    }else {RrO_h=NULL}
    RrO_h
  })
  return(list(RrE,RrO))
}

# from the output of the constraints in 'parse_hypothesis' create lists for the equality and order matrices
# different format parse_hyp object
make_RrList2 <- function(parse_hyp2){
  numhyp <- length(parse_hyp2$original_hypothesis)
  qE <- parse_hyp2$n_constraints[(0:(numhyp-1))*2+1]
  qO <- parse_hyp2$n_constraints[(1:numhyp)*2]
  RrE <- lapply(1:numhyp,function(h){
    startcon <- sum(qE[1:h]+qO[1:h])-qE[h]-qO[h]
    if(qE[h]==1){
      RrE_h <- t(as.matrix(parse_hyp2$hyp_mat[startcon+1:qE[h],]))
    }else if(qE[h]>1){
      RrE_h <- parse_hyp2$hyp_mat[startcon+1:qE[h],]
    }else {RrE_h=NULL}
    RrE_h
  })
  RrO <- lapply(1:numhyp,function(h){
    startcon <- sum(qE[1:h]+qO[1:h])-qE[h]-qO[h]
    if(qO[h]==1){
      RrO_h <- t(as.matrix(parse_hyp2$hyp_mat[startcon+qE[h]+1:qO[h],]))
    }else if(qO[h]>1){
      RrO_h <- parse_hyp2$hyp_mat[startcon+qE[h]+1:qO[h],]
    }else {RrO_h=NULL}
    RrO_h
  })
  return(list(RrE,RrO))
}

#for checking whether constraints are conflicting replace interval constraints by equality constraints
interval_RrStack <- function(RrStack){
  q1 <- nrow(RrStack)
  q2 <- ncol(RrStack)
  RrStack_out <- RrStack
  if(q1 > 1){
    row1 <- 1
    while(row1 < q1){
      for(row2 in (row1+1):q1){
        #        print(row2)
        if(sum(abs(RrStack_out[row1,-q2] + RrStack_out[row2,-q2]))==0){ # && RrStack_out[row1,q2]!=RrStack_out[row2,q2] ){
          #together row1 and row2 imply an interval constraint
          whichcol <- abs(RrStack_out[row1,-q2])!=0
          whichcol1 <- which(whichcol)
          if(sum(whichcol)==1){
            welkpos <- ifelse(RrStack_out[row1,c(whichcol,F)]>0,row1,row2)
            welkneg <- ifelse(RrStack_out[row1,c(whichcol,F)]<0,row1,row2)
            lb <- RrStack_out[welkpos,q2]
            ub <- -RrStack_out[welkneg,q2]
            RrStack_out[row1,] <- RrStack_out[welkpos,]
            RrStack_out[row1,q2] <- (ub+lb)/2
            RrStack_out <- RrStack_out[-row2,]
            q1 <- q1 - 1
          }else{
            RrStack_out[row1,q2] <- 0
            RrStack_out <- RrStack_out[-row2,]
            q1 <- q1 - 1
          }
          break
        }
      }
      row1 <- row1 + 1
    }
  }
  if(is.matrix(RrStack_out)==F){
    RrStack_out <- t(RrStack_out)
  }
  return(RrStack_out)
}

params_in_hyp <- function(hyp){
  params_in_hyp <- trimws(unique(strsplit(hyp, split = "[ =<>,\\(\\);&\\*+-]+", perl = TRUE)[[1]]))
  params_in_hyp <- params_in_hyp[!sapply(params_in_hyp, grepl, pattern = "^[0-9]*\\.?[0-9]+$")]
  params_in_hyp[grepl("^[a-zA-Z]", params_in_hyp)]
}
