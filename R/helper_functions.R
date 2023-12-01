check_vcov <- function(x){
  if (!isTRUE(all.equal(x, t(x))) || any(diag(x) < 0)){
    saveRDS(x, "c:/git_repositories/BFpack/erordump.RData")
    stop(sQuote("sigma"), " is not a covariance matrix")
  }
}

process.prior.hyp.explo <- function(prior_hyp_explo, model){

  if(sum(class(model)=="bartlett_htest")>0){
    if(is.null(prior_hyp_explo)){
      prior_hyp_explo <- rep(.5,2)
    }else{
      if(length(prior_hyp_explo)!=2){
        stop("For an object of class 'bartlett_htest', the argument 'prior_hyp_explo' should be a vector
             of length 2 of the prior probabilities for the hypotheses of homogeneity of variances and an
             unconstrained alternative. Or use the default ('NULL') which implies equal prior probabilities.")
      }
    }
  }else{
    if(is.null(prior_hyp_explo)){ # then equal prior probabilities
      if(sum(class(model)=="aov")>0){
        prior_hyp_explo <- list(rep(1,3),rep(1,2),rep(1,2))
      }else{
        prior_hyp_explo <- list(rep(1,3))
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
            prior_hyp_explo <- list(rep(1,3),rep(1,2),rep(1,2))
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
                prior_hyp_explo[[1]] <- rep(1,3)
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
                prior_hyp_explo[[1]] <- rep(1,3)
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
              prior_hyp_explo <- list(rep(1,3))
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


