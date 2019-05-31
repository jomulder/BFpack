
#' @method print BF
#' @export
print.BF <- function(x, #stats = c("Fit_eq", "Com_eq", "Fit_in", "Com_in", "Fit", "Com", "BF", "PMPa", "PMPb"),
                       digits = 3,
                       na.print = "", ...){

  cat("Results of Bayesian exploratory hypothesis tests for an object of class ", class(x$model)[1], ":\n\n", sep = "")
  cat("Coefficients:")
  cat("\n")
  print(round(x$PHP_exploratory,digits))

  cat("\n")

  if(!is.null(x$BFtu_confirmatory)){
    cat("Results of Bayesian confirmatory hypothesis tests for an object of class ", class(x$model)[1], ":\n\n", sep = "")
    cat("Hypotheses:")
    cat("\n")

    PHPmatrix <- as.matrix(round(x$PHP_confirmatory,digits))
    colnames(PHPmatrix) <- "Pr(hypothesis)"
    print(PHPmatrix)

    cat("\n")
    cat("Evidence matrix:")
    cat("\n")

    print(round(x$BFmatrix_confirmatory,digits))
  }
  # cat("\nHypotheses:\n ", paste(rownames(dat)[-nrow(dat)], ": ", x$hypotheses, sep = "", collapse = "\n  "))
  #
  # if(!is.null(x[["warnings"]])){
  #   warning("Bain analysis returned the following warnings:\n  ", paste(1:length(x$warnings), ". ", x$warnings, sep = "", collapse = "\n  "))
  # }
}




