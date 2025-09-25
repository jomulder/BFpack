
#' @method print BF
#' @export
print.BF <- function(x,
                     digits = 3,
                     na.print = "", ...){


  cat("Call:")
  cat("\n")
  print(x$call)

  cat("\n")

  digits <- 3

  if(is.null(x$BFtu_confirmatory)){

    cat("Bayesian hypothesis test","\n", sep = "")
    cat("Type: exploratory","\n", sep = "")
    cat("Object: ",class(x$model)[1],"\n", sep = "")
    cat("Parameter: ",x$parameter,"\n", sep = "")
    cat("Method: ",x$bayesfactor,"\n\n", sep = "")

    cat("Posterior probabilities of the hypotheses","\n", sep = "")
    cat("\n")
    cat("Model parameters:")
    cat("\n")
    print(round(x$PHP_exploratory,digits))

    cat("\n")

    if(sum(class(x$model)=="aov")>0){
      if(!is.null(x$PHP_main)){
        cat("main effects:")
        cat("\n")
        print(round(x$PHP_main,digits))

        cat("\n")
      }
      if(!is.null(x$PHP_interaction)){
        cat("interaction effects:")
        cat("\n")
        print(round(x$PHP_interaction,digits))

        cat("\n")
      }
    }
    if(class(x$model)[1]=="rma.uni"){
      cat("Posterior support for random effects","\n", sep = "")
      cat("\n")
      tau2out <- matrix(x$estimates[2,5])
      row.names(tau2out) <- "tau2 (marema)"
      colnames(tau2out) <- "Pr(>0)"
      print(round(tau2out,digits))
      cat("\n")
    }
    if(class(x$model)[1]=="cor_test"){
      cat("Bartlett's test of sphericity (all correlations are zero):","\n", sep = "")
      cat("\n")
      if(nrow(x$PHP_Bartlett)>1){
        print(round(x$PHP_Bartlett,digits))
      }else{
        print(round(x$PHP_Bartlett[1,],digits))
      }
      cat("\n")
    }

  }else{

    cat("Bayesian hypothesis test","\n", sep = "")
    cat("Type: confirmatory","\n", sep = "")
    cat("Object: ",class(x$model)[1],"\n", sep = "")
    cat("Parameter: ",x$parameter,"\n", sep = "")
    cat("Method: ",x$bayesfactor,"\n\n", sep = "")

    cat("Posterior probabilities of the hypotheses")
    cat("\n")

    PHPmatrix <- as.matrix(round(x$PHP_confirmatory,digits))
    colnames(PHPmatrix) <- "Pr(hypothesis|data)"
    hypnumbers <- unlist(lapply(1:nrow(PHPmatrix),function(r){
      paste0("H",as.character(r))
    }))
    row.names(PHPmatrix) <- hypnumbers
    print(PHPmatrix)

    cat("\n")
    if(x$log==FALSE){
      cat("Bayes factors (rows versus columns):")
    }else{
      cat("Log Bayes factors (rows versus columns):")
    }
    cat("\n")

    BFmat <- round(x$BFmatrix_confirmatory,digits)
    row.names(BFmat) <- colnames(BFmat) <- hypnumbers
    print(BFmat)

    cat("\n")
    cat("Hypotheses:")
    cat("\n")

    for(h in 1:length(x$hypotheses)){
      cat(paste0(hypnumbers[h],": ",x$hypotheses[h]))
      cat("\n")
    }

  }

}




