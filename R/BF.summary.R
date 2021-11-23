
#' @method summary BF
#' @export
summary.BF <- function(object, ...){

  cat("Call:")
  cat("\n")
  print(object$call)

  cat("\n")
  digits <- 3

  cat("Bayesian hypothesis test","\n", sep = "")
  cat("Type: exploratory","\n", sep = "")
  cat("Object: ",class(object$model)[1],"\n", sep = "")
  cat("Parameter: ",object$parameter,"\n", sep = "")
  cat("Method: ",object$bayesfactor,"\n\n", sep = "")

  cat("Posterior probabilities:","\n", sep = "")
  print(round(object$PHP_exploratory,digits))

  cat("\n")

  if(!is.null(object$BFtu_confirmatory)){

    cat("Bayesian hypothesis test","\n", sep = "")
    cat("Type: confirmatory","\n", sep = "")
    cat("Object: ",class(object$model)[1],"\n", sep = "")
    cat("Parameter: ",object$parameter,"\n", sep = "")
    cat("Method: ",object$bayesfactor,"\n\n", sep = "")

    cat("Posterior probabilities:")
    cat("\n")

    PHPmatrix <- as.matrix(round(object$PHP_confirmatory,digits))
    colnames(PHPmatrix) <- "Pr(hypothesis|data)"
    hypnumbers <- unlist(lapply(1:nrow(PHPmatrix),function(r){
      paste0("H",as.character(r))
    }))
    row.names(PHPmatrix) <- hypnumbers
    print(PHPmatrix)

    cat("\n")
    cat("Evidence matrix (Bayes factors):")
    cat("\n")

    BFmat <- round(object$BFmatrix_confirmatory,digits)
    row.names(BFmat) <- colnames(BFmat) <- hypnumbers
    print(BFmat)

    cat("\n")
    cat("Specification table:")
    cat("\n")

    BFtable <- round(object$BFtable_confirmatory,digits)
    row.names(BFtable) <- hypnumbers
    print(BFtable)

    cat("\n")
    cat("Hypotheses:")
    cat("\n")

    for(h in 1:length(object$hypotheses)){
      cat(paste0(hypnumbers[h],": ",object$hypotheses[h]))
      cat("\n")
    }
  }
}




