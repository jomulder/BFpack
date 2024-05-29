#' @method print cor_test
#' @export
print.cor_test <- function(x,
                           digits = 3, ...){

  estimates <- x$correstimates
  names <- x$corrnames
  groups <- length(names)
  P <- nrow(names[[1]])
  numcorr <- P*(P-1)/2
  countg = 0
  corrlist <- lapply(1:groups,function(g){
    lapply(1:3,function(b){
      matje <- matrix(NA,P,P)
      row.names(matje) <- colnames(matje) <- x$variables[[1]]
      matje[lower.tri(diag(P))] <- estimates[numcorr*(g-1)+1:numcorr,1+b]
      matje
    })
  })

  cat("\n")
  cat("Unconstrained Bayesian estimates","\n", sep = "")
  cat("\n")

  if(groups > 1){
    for(g in 1:groups){
      cat(paste0("Group g",as.character(g),":"),"\n", sep = "")
      cat("\n")
      cat("correlation types","\n")
      print(x$cor.type[[g]],na.print="",quote=FALSE)
      cat("\n")
      cat("Posterior 2.5% lower bounds:","\n", sep = "")
      print(round(corrlist[[g]][[2]],digits), na.print = "")
      cat("\n")
      cat("Posterior median:","\n", sep = "")
      print(round(corrlist[[g]][[1]],digits), na.print = "")
      cat("\n")
      cat("Posterior 97.5% upper bounds:","\n", sep = "")
      print(round(corrlist[[g]][[3]],digits), na.print = "")
      cat("\n")
    }
  }else{
    cat("correlation types","\n")
    print(x$cor.type[[1]],na.print="",quote=FALSE)
    cat("\n")
    cat("Posterior 2.5% lower bounds:","\n", sep = "")
    print(round(corrlist[[1]][[2]],digits), na.print = "")
    cat("\n")
    cat("Posterior median:","\n", sep = "")
    print(round(corrlist[[1]][[1]],digits), na.print = "")
    cat("\n")
    cat("Posterior 97.5% upper bounds:","\n", sep = "")
    print(round(corrlist[[1]][[3]],digits), na.print = "")
    cat("\n")
  }

}


#' @method summary cor_test
#' @export
summary.cor_test <- function(x,
                           digits = 3, ...){

  cor.df <- round(as.data.frame(x$correstimates),digits)
  cor.df$cor.type <- unlist(lapply(1:length(x$cor.type),function(g){
    x$cor.type[[g]][lower.tri(x$cor.type[[g]])]
  }))

  cor.df

}
