
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

    cat("Posterior probabilities:","\n", sep = "")
    print(round(x$PHP_exploratory,digits))

    cat("\n")

  }else{

    cat("Bayesian hypothesis test","\n", sep = "")
    cat("Type: confirmatory","\n", sep = "")
    cat("Object: ",class(x$model)[1],"\n", sep = "")
    cat("Parameter: ",x$parameter,"\n", sep = "")
    cat("Method: ",x$bayesfactor,"\n\n", sep = "")

    cat("Posterior probabilities:")
    cat("\n")

    PHPmatrix <- as.matrix(round(x$PHP_confirmatory,digits))
    colnames(PHPmatrix) <- "Pr(hypothesis|data)"
    hypnumbers <- unlist(lapply(1:nrow(PHPmatrix),function(r){
      paste0("H",as.character(r))
    }))
    row.names(PHPmatrix) <- hypnumbers
    print(PHPmatrix)

    cat("\n")
    cat("Evidence matrix:")
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




