
#' @method summary BF
#' @export
summary.BF <- function(object, ...){

  cat("Call:")
  cat("\n")
  print(object$call)

  cat("\n")
  digits <- 3

  cat("BFpack: Exploratory Bayes factor tests for an object of class ", class(object$model)[1], ":\n\n", sep = "")
  # cat("Coefficients:")
  # cat("\n")
  print(round(object$PHP_exploratory,digits))

  cat("\n")

  if(class(object$model)[1]=="aov"){
    if(!is.null(object$BFtu_main)){
      cat("Main effects:")
      cat("\n")
      print(round(object$PHP_main,digits))
    }
    cat("\n")

    if(!is.null(object$BFtu_interaction)){
      cat("Interaction effects:")
      cat("\n")
      print(round(object$PHP_interaction,digits))
    }
    cat("\n")
  }

  if(!is.null(object$BFtu_confirmatory)){
    cat("BFpack: Confirmatory Bayes factor tests for an object of class ", class(object$model)[1], ":\n\n", sep = "")
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
    cat("Evidence matrix:")
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

    # hyps <- as.matrix(x$hypotheses)
    # row.names(hyps) <- hypnumbers
    # colnames(hyps) <- NULL
    # print(hyps)
    for(h in 1:length(object$hypotheses)){
      # hyp <- x$hypotheses[h]
      # names(hyp) <- hypnumbers[h]
      # print(hyp)
      cat(paste0(hypnumbers[h],": ",object$hypotheses[h]))
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




