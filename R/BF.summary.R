
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

  cat("Posterior probabilities of the hypotheses","\n", sep = "")
  cat("\n")
  cat("Model parameters:")
  cat("\n")
  print(round(object$PHP_exploratory,digits))

  cat("\n")

  if(sum(class(object$model)=="aov")>0){
    if(!is.null(object$PHP_main)){
      cat("main effects:")
      cat("\n")
      print(round(object$PHP_main,digits))

      cat("\n")
    }
    if(!is.null(object$PHP_interaction)){
      cat("interaction effects:")
      cat("\n")
      print(round(object$PHP_interaction,digits))

      cat("\n")
    }
  }
  if(class(object$model)[1]=="rma.uni"){
    cat("Posterior support for random effects","\n", sep = "")
    cat("\n")
    tau2out <- matrix(object$estimates[2,5])
    row.names(tau2out) <- "tau2 (marema)"
    colnames(tau2out) <- "Pr(>0)"
    print(round(tau2out,digits))
    cat("\n")
  }
  if(class(object$model)[1]=="cor_test"){
    cat("Bartlett's test of sphericity (all correlations are zero):","\n", sep = "")
    cat("\n")
    if(nrow(object$PHP_Bartlett)>1){
      print(round(object$PHP_Bartlett,digits))
    }else{
      print(round(object$PHP_Bartlett[1,],digits))
    }
    cat("\n")
  }

  if(!is.null(object$BFtu_confirmatory)){

    cat("Bayesian hypothesis test","\n", sep = "")
    cat("Type: confirmatory","\n", sep = "")
    cat("Object: ",class(object$model)[1],"\n", sep = "")
    cat("Parameter: ",object$parameter,"\n", sep = "")
    cat("Method: ",object$bayesfactor,"\n\n", sep = "")

    cat("Posterior probabilities of the hypotheses")
    cat("\n")

    PHPmatrix <- as.matrix(round(object$PHP_confirmatory,digits))
    colnames(PHPmatrix) <- "Pr(hypothesis|data)"
    hypnumbers <- unlist(lapply(1:nrow(PHPmatrix),function(r){
      paste0("H",as.character(r))
    }))
    row.names(PHPmatrix) <- hypnumbers
    print(PHPmatrix)

    cat("\n")
    if(object$log==FALSE){
      cat("Bayes factors (rows versus columns):")
    }else{
      cat("Log Bayes factors (rows versus columns):")
    }
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




