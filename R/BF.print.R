
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
#    cat("BFpack: Exploratory Bayes factor tests for an object of class ", class(x$model)[1], ":\n\n", sep = "")

    cat("Bayesian hypothesis test","\n", sep = "")
    cat("Type: Exploratory","\n", sep = "")
    cat("Object: ",class(x$model)[1],"\n", sep = "")
    cat("Parameter: ",x$parameter,"\n", sep = "")
    cat("Method: ",x$bayesfactor,"\n\n", sep = "")
    # cat("\n")
    cat("Posterior probabilities:","\n", sep = "")
    print(round(x$PHP_exploratory,digits))

    cat("\n")

    if(class(x$model)[1]=="aov"){
      if(!is.null(x$BFtu_main)){
        cat("Main effects:")
        cat("\n")
        print(round(x$PHP_main,digits))
      }
      cat("\n")

      if(!is.null(x$BFtu_interaction)){
        cat("Interaction effects:")
        cat("\n")
        print(round(x$PHP_interaction,digits))
      }
      cat("\n")
    }
  }else{

    cat("Bayesian hypothesis test","\n", sep = "")
    cat("Type: Confirmatory","\n", sep = "")
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

    # cat("\n")
    # cat("Evidence matrix:")
    # cat("\n")
    #
    # BFmat <- round(x$BFmatrix_confirmatory,digits)
    # row.names(BFmat) <- colnames(BFmat) <- hypnumbers
    # print(BFmat)

    cat("\n")
    cat("Hypotheses:")
    cat("\n")

    # hyps <- as.matrix(x$hypotheses)
    # row.names(hyps) <- hypnumbers
    # colnames(hyps) <- NULL
    # print(hyps)
    for(h in 1:length(x$hypotheses)){
      # hyp <- x$hypotheses[h]
      # names(hyp) <- hypnumbers[h]
      # print(hyp)
      cat(paste0(hypnumbers[h],": ",x$hypotheses[h]))
      cat("\n")
    }
    # hypfull <- unlist(lapply(1:length(x$hypotheses),function(r){
    #   paste0("H",as.character(r),": ",x$hypotheses[r])
    # }))

  }
  # cat("\nHypotheses:\n ", paste(rownames(dat)[-nrow(dat)], ": ", x$hypotheses, sep = "", collapse = "\n  "))
  #
  # if(!is.null(x[["warnings"]])){
  #   warning("Bain analysis returned the following warnings:\n  ", paste(1:length(x$warnings), ". ", x$warnings, sep = "", collapse = "\n  "))
  # }
}




