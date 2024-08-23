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
  if(x$prior.cor=="joint.unif"){
    cat("Unconstrained Bayesian estimation (joint uniform prior)","\n", sep = "")
  }else{
    cat("Unconstrained Bayesian estimation (marginally uniform prior)","\n", sep = "")
  }
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
summary.cor_test <- function(object, digits = 3, ...){

  cor.df <- round(as.data.frame(object$correstimates),digits)
  cor.df$cor.type <- unlist(lapply(1:length(object$cor.type),function(g){
    object$cor.type[[g]][lower.tri(object$cor.type[[g]])]
  }))

  cor.df

}


#' @importFrom coda mcmc
#' @method plot cor_test
#' @export
plot.cor_test <- function(x, ...){

  numgroups <- length(x$corrdraws)
  P <- dim(x$corrdraws[[1]])[2]
  numcor <- P*(P-1)/2
  numcor.total <- numcor * numgroups

  cor.draws.matrix <- mcmc(data=do.call(cbind,lapply(1:numgroups,function(g){
    do.call(cbind,lapply(1:(P-1),function(p){
      draws_g_p <- as.matrix(x$corrdraws[[g]][,(p+1):P,p])
      colnames(draws_g_p) <- x$corrnames[[g]][(p+1):P,p]
      draws_g_p
    }))
  })),start=1,end=length(x$corrdraws[[1]][,1,1]),thin=1)

  plot(cor.draws.matrix, ...)

}



