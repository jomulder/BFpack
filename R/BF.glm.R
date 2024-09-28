#BF method for glm classes


#' @method BF glm
#' @export
BF.glm <- function(x,
                   hypothesis = NULL,
                   prior.hyp.explo = NULL,
                   prior.hyp.conf = NULL,
                   prior.hyp = NULL,
                   complement = TRUE,
                   log = FALSE,
                   cov.prob = .95,
                   ...){

  if(as.character(x$family)[1]=="gaussian"){

    # then use BF.lm
    class(x) <- "lm"
    Args <- as.list(match.call()[-1])
    Args$x <- x
    Args$hypothesis <- hypothesis
    Args$prior.hyp <- prior.hyp
    Args$prior.hyp.conf <- prior.hyp.conf
    Args$prior.hyp.explo <- prior.hyp.explo
    Args$complement <- complement
    Args$log <- log
    Args$cov.prob <- cov.prob
    out <- do.call(BF, Args)
    out$model <- x
    out$call <- match.call()

  }else{

    Args <- as.list(match.call()[-1])
    get_est <- get_estimates(x)
    Args$x <- get_est$estimate
    Args$Sigma <- get_est$Sigma[[1]]
    Args$n <- nobs(x)
    Args$hypothesis <- hypothesis
    Args$prior.hyp <- prior.hyp
    Args$prior.hyp.conf <- prior.hyp.conf
    Args$prior.hyp.explo <- prior.hyp.explo
    Args$complement <- complement
    Args$log <- log
    out <- do.call(BF, Args)
    out$model <- x
    out$call <- match.call()

  }

  out
}


