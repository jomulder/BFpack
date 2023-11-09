#' @title Bartlett Test of Homogeneity of Variances
#' @description Performs Bartlett's test of the null that the
#' variances in each of the groups (samples) are the same.
#'
#'@details \code{x} must be a numeric data vector, and \code{g}
#'must be a vector or factor object of the same length as \code{x}
#'giving the group for the corresponding elements of \code{x}.
#'
#'@section Bain t_test:
#'In order to allow users to enjoy the functionality of bain with the familiar
#'stats-function \code{bartlett.test}, we have had to make minor changes to the
#'function \code{bartlett.test.default}. All rights to, and credit for, the
#'function \code{bartlett.test.default}
#'belong to the R Core Team, as indicated in the original license below.
#'We make no claims to copyright and incur no liability with regard to the
#'changes implemented in \code{bartlett_test}.
#'
#'This the original copyright notice by the R core team:
#'File src/library/stats/R/bartlett_test.R
#'Part of the R package, https://www.R-project.org
#'
#'Copyright (C) 1995-2015 The R Core Team
#'
#' This program is free software; you can redistribute it and/or modify
#' it under the terms of the GNU General Public License as published by
#' the Free Software Foundation; either version 2 of the License, or
#' (at your option) any later version.
#'
#' This program is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU General Public License for more details.
#'
#' A copy of the GNU General Public License is available at
#' https://www.R-project.org/Licenses/
#'
#'@aliases bartlett_test bartlett_test.default
#'@param x a numeric vector of data values, or a list of
#'numeric data vectors representing the respective samples,
#'or fitted linear model objects (inheriting from class "lm").
#'@param g a vector or factor object giving the group for
#'the corresponding elements of x. Ignored if x is a list.
#'@param ... further arguments to be passed to or from methods.
#'
#'@return A list with class \code{"bartlett_htest"} containing the following
#'components: \item{statistic}{Bartlett's K-squared test statistic.}
#'\item{parameter}{the degrees of freedom of the approximate chi-squared
#'distribution of the test statistic.}
#'\item{p.value}{the p-value of the test.} \item{conf.int}{a confidence
#'interval for the mean appropriate to the specified alternative hypothesis.}
#'\item{method}{the character string "Bartlett test of homogeneity of variances".}
#'\item{data.name}{a character string giving the names of the data.}
#'\item{vars}{the sample variances across groups (samples).}
#'\item{n}{the number of observations per group (sample)}
#'
#'@references Bartlett, M. S. (1937). Properties of sufficiency
#'and statistical tests. Proceedings of the Royal Society of
#'London Series A 160, 268–282. DOI: 10.1098/rspa.1937.0109.
#'
#'@examples
#'require(graphics)
#'
#'plot(count ~ spray, data = InsectSprays)
#'bartlett_test(InsectSprays$count, InsectSprays$spray)
#'
#' @rdname bartlett_test
#' @export
bartlett_test <- function(x, g, ...) UseMethod("bartlett_test", x)


#' @importFrom stats bartlett.test
#' @method bartlett_test default
#' @rdname bartlett_test
#' @export
bartlett_test.default <- function(x, g, ...){
  cl <- match.call()
  cl[[1]] <- as.name("bartlett.test")
  bart <- eval.parent(cl)
  vars <- tapply(x, g, var, na.rm = TRUE)
  n <- table(g)
  names(vars) <- names(n)
  class(bart) <- c("bartlett_htest", class(bart))
  bart$vars <- vars
  bart$n <- n
  bart
}

# require(stats)
# exists("bartlett.test.default") # false
# getS3method("bartlett.test", "default")

#' @importFrom stats rgamma
#' @importFrom stats rchisq
#' @method BF bartlett_htest
#' @export
BF.bartlett_htest <- function(x,
                           hypothesis = NULL,
                           prior.hyp = NULL,
                           complement = TRUE,
                           log = FALSE,
                           ...) {
  get_est <- get_estimates(x)
  nsim <- 1e5
  s2 <- get_est$estimate
  n <- c(x$n)
  b <- 2/n
  J <- length(n)
  names_coef <- names(get_est$estimate)

  logIN <- log

  # exploratory BF for equality of variances:
  logmx0 <- - 1 / 2 * sum((1 - b) * n) * log(pi) + 1 / 2 * log(prod(b)) +
    lgamma((sum(n) - J) / 2) - lgamma((sum(b * n) - J) / 2) -
    1 / 2 * (sum(n) - J) * log(sum((n - 1) * s2)) +
    1 / 2 * (sum(b * n) - J) * log(sum(b * (n - 1) * s2))
  logmxu <- - 1 / 2 * sum((1 - b) * n) * log(pi) + 1 / 2 * log(prod(b)) +
    sum(lgamma((n - 1) / 2) - lgamma((b * n - 1) / 2) -
          1 / 2 * (n - 1) * log((n - 1) * s2) +
          1 / 2 * (b * n - 1) * log(b * (n - 1) * s2))
  BF0u <- logmx0 - logmxu

  BFtu_exploratory <- c(BF0u,log(1))
  names(BFtu_exploratory) <- c("homogeneity of variances","no homogeneity of variances")
  PHP_exploratory <- exp(BFtu_exploratory - max(BFtu_exploratory)) /
    sum(exp(BFtu_exploratory - max(BFtu_exploratory)))

  if (!is.null(hypothesis)){
    parse_hyp <- parse_hypothesis(names_coef, hypothesis)
    parse_hyp$hyp_mat <- do.call(rbind, parse_hyp$hyp_mat)
    RrList <- make_RrList2(parse_hyp)
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]
  }

  if (is.null(hypothesis)) {
    BFmatrix_confirmatory <- PHP_confirmatory <- BFtu_confirmatory <- relfit <-
      relcomp <- hypotheses <- BFtable <- priorprobs <- NULL
  } else if (all(unlist(lapply(append(RrE, RrO), is.null)))) {
    BFmatrix_confirmatory <- PHP_confirmatory <- BFtu_confirmatory <- relfit <-
      relcomp <- hypotheses <- BFtable <- priorprobs <- NULL
  } else { # execute confirmatory Bayes factor test based on hypothesis input
    # check if hypotheses are admissible:
    RrCheck <- do.call(rbind, append(RrE, RrO))
    RrCheck_count <- t(apply(RrCheck[, -ncol(RrCheck), drop = FALSE], 1,
                             function(x) {sapply(list(-1, 1), function(y) {sum(y == x)})}))
    if (any(RrCheck_count != 1) || any(RrCheck[, ncol(RrCheck)] != 0)) {
      stop(paste0("The hypotheses contain inadmissible constraints."))
    }

    Th <- length(RrE)
    logmx <- relfit <- relcomp <- logmxE <- rep(NA, times = Th)
    names(logmx) <- names(relfit) <- names(relcomp) <- names(logmxE) <-
      parse_hyp$original_hypothesis

    for (h in 1:Th) {

      if (is.null(RrE[[h]])) {
        unique_vars <- as.list(1:J)
      } else {
        RrEh <- RrE[[h]][, -ncol(RrE[[h]])]
        if (!is.matrix(RrEh)) {
          RrEh <- t(as.matrix(RrEh))
        }
        RrEh_pos <- t(apply(RrEh, 1, function(x) which(!(x == 0))))
        unique_vars <- list()
        rows <- 1:nrow(RrEh_pos)
        while (length(rows) > 0) {
          equal_vars <- RrEh_pos[min(rows), ]
          row_check <- min(rows)
          for (i in setdiff(rows, row_check)) {
            if (any(equal_vars %in% RrEh_pos[i, ])) {
              equal_vars <- unique(c(equal_vars, RrEh_pos[i, ]))
              row_check <- c(row_check, i)
            }
          }
          unique_vars <- c(unique_vars, list(equal_vars))
          rows <- setdiff(rows, row_check)
        }
        unique_vars <- c(unique_vars, setdiff(1:J, unlist(unique_vars)))
      }
      K <- length(unique_vars)
      Jk <- sapply(unique_vars, length)
      s2list <- lapply(unique_vars, function(x) s2[x])
      nlist  <- lapply(unique_vars, function(x) n[x])
      blist  <- lapply(unique_vars, function(x) b[x])
      df <- dfb <- SS <- SSb <- rep(NA, times = K)
      for (i in 1:K) {
        df[i]  <- sum(nlist[[i]]) - Jk[i]
        dfb[i] <- sum(blist[[i]] * nlist[[i]]) - Jk[i]
        SS[i]  <- sum((nlist[[i]] - 1) * s2list[[i]])
        SSb[i] <- sum(blist[[i]] * (nlist[[i]] - 1) * s2list[[i]])
      }
      logmxE[h] <- - 1 / 2 * sum((1 - unlist(blist)) * unlist(nlist)) * log(pi) +
        1 / 2 * log(prod(unlist(blist))) + sum(lgamma(df / 2) - lgamma(dfb / 2) -
                                                 1 / 2 * df * log(SS) + 1 / 2 * dfb * log(SSb))
      if (is.null(RrO[[h]])) {
        logmx[h] <- logmxE[h]
      } else {
        RrOh <- RrO[[h]][, -ncol(RrO[[h]])]
        if (!is.matrix(RrOh)) {
          RrOh <- t(as.matrix(RrOh))
        }
        RrOh_pos <- t(apply(RrOh, 1, function(x) c(which(x == -1), which(x == 1))))
        unique_vars_order <- cbind(
          apply(as.matrix(RrOh_pos[, 1]), 1, function(x) {
            which(unlist(lapply(unique_vars, function(y) {x %in% y})))}),
          apply(as.matrix(RrOh_pos[, 2]), 1, function(x) {
            which(unlist(lapply(unique_vars, function(y) {x %in% y})))})
        )
        post_samp <- prior_samp <- matrix(NA, nrow = nsim, ncol = K)
        indi_post <- indi_prior <- rep(1, times = nsim)
        for (i in unique(c(unique_vars_order))) {
          post_samp[, i] <- SS[i] / rchisq(nsim, df = df[i])
          prior_samp[, i] <- dfb[i] / rchisq(nsim, df = dfb[i])
        }
        for (i in 1:nrow(unique_vars_order)) {
          indi_post <- indi_post * (post_samp[, unique_vars_order[i, 1]] <
                                      post_samp[, unique_vars_order[i, 2]])
          indi_prior <- indi_prior * (prior_samp[, unique_vars_order[i, 1]] <
                                        prior_samp[, unique_vars_order[i, 2]])
        }
        relfit[h] <- sum(indi_post) / nsim
        relcomp[h] <- sum(indi_prior) / nsim
        logmx[h] <- log(relfit[h] / relcomp[h]) + logmxE[h]
      }
    }
    if(complement==TRUE){
      #compute marginal likelihood for complement hypothesis
      relfit <- inversegamma_prob_Hc(shape1=(n-1)/2,scale1=s2*(n-1)/(2*n),relmeas=relfit,RrE1=RrE,RrO1=RrO)
      relcomp <- inversegamma_prob_Hc(shape1=rep(.5,length(n)),scale1=rep(.5,length(n)),relmeas=relcomp,RrE1=RrE,RrO1=RrO)
      if(length(relfit)>Th){
        logmxE <- c(logmxE,logmxu)
        logmx <- c(logmx,logmxu + log(relfit[Th+1]/relcomp[Th+1]))
        names(logmx)[Th+1] <- "complement"
      }
    }
    hypotheses <- names(logmx)
    BFtu_confirmatory <- (logmx - logmxu)
    BFmatrix_confirmatory <- BFtu_confirmatory %*% t(rep(1,length(BFtu_confirmatory))) -
      rep(1,length(BFtu_confirmatory)) %*% t(BFtu_confirmatory)
    diag(BFmatrix_confirmatory) <- log(1)
    names(BFtu_confirmatory) <- row.names(BFmatrix_confirmatory) <-
      colnames(BFmatrix_confirmatory) <- hypotheses

    if(is.null(prior.hyp)){
      priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
    }else{
      if(!is.numeric(prior.hyp) || length(prior.hyp)!=length(BFtu_confirmatory)){
        warning(paste0("Argument 'prior.hyp' should be numeric and of length ",as.character(length(BFtu_confirmatory)),". Equal prior probabilities are used."))
        priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
      }else{
        priorprobs <- prior.hyp
      }
    }

    PHP_confirmatory <- exp(BFtu_confirmatory - max(BFtu_confirmatory)) * priorprobs /
      sum(exp(BFtu_confirmatory - max(BFtu_confirmatory)) * priorprobs)
    relcomp[which(is.na(relcomp))] <- 1
    relfit[which(is.na(relfit))] <- 1
    BF_E <- exp(logmxE - logmxu)
    BFtable <- cbind(rep(NA,length(relfit)),relcomp,rep(NA,length(relfit)),relfit,BF_E,
                     relfit/relcomp,BF_E*relfit/relcomp,PHP_confirmatory)
    row.names(BFtable) <- names(PHP_confirmatory)
    colnames(BFtable) <- c("complex=","complex>","fit=","fit>","BF=","BF>","BF","PHP")

    if(logIN == FALSE){
      BFtu_exploratory <- exp(BFtu_exploratory)
      BFtu_confirmatory <- exp(BFtu_confirmatory)
      BFmatrix_confirmatory <- exp(BFmatrix_confirmatory)
    }
  }

  BFlm_out <- list(
    BFtu_exploratory=BFtu_exploratory,
    PHP_exploratory=PHP_exploratory,
    BFtu_confirmatory=BFtu_confirmatory,
    PHP_confirmatory=PHP_confirmatory,
    BFmatrix_confirmatory=BFmatrix_confirmatory,
    BFtable_confirmatory=BFtable,
    prior.hyp=priorprobs,
    hypotheses=hypotheses,
    estimates=s2,
    model=x,
    bayesfactor="generalized adjusted fractional Bayes factors",
    parameter="group variances",
    log = logIN,
    call=match.call())

  class(BFlm_out) <- "BF"

  return(BFlm_out)
}


# The function computes the probability of the complement hypothesis
inversegamma_prob_Hc <- function(shape1,scale1,relmeas,RrE1,RrO1,samsize1=1e5){

  numhyp <- length(RrE1)
  whichO <- unlist(lapply(1:numhyp,function(h){is.null(RrE1[[h]])}))
  numO <- sum(whichO)
  numpara <- length(shape1)

  if(numO==length(RrE1)){ # Then the complement is equivalent to the unconstrained hypothesis.
    relmeas <- c(relmeas,1)
    names(relmeas)[numhyp+1] <- "complement"
  }else{ # So there is at least one hypothesis with only order constraints
    if(numO==1){ # There is one hypothesis with only order constraints. Hc is complement of this hypothesis.
      relmeas <- c(relmeas,1-relmeas[whichO])
      names(relmeas)[numhyp+1] <- "complement"
    }else{ # So more than one hypothesis with only order constraints
      randomDraws <- rmvnorm(samsize1,mean=rep(0,numpara),sigma=diag(numpara))
      #get draws that satisfy the constraints of the separate order constrained hypotheses
      checksOC <- lapply(which(whichO),function(h){
        Rorder <- as.matrix(RrO1[[h]][,-(1+numpara)])
        if(ncol(Rorder)==1){
          Rorder <- t(Rorder)
        }
        rorder <- as.matrix(RrO1[[h]][,1+numpara])
        apply(randomDraws%*%t(Rorder) > rep(1,samsize1)%*%t(rorder),1,prod)
      })
      checkOCplus <- Reduce("+",checksOC)

      if(sum(checkOCplus > 0) < samsize1){ #then the joint order constrained hypotheses do not completely cover the parameter space.
        if(sum(checkOCplus>1)==0){ # then order constrained spaces are nonoverlapping
          relmeas <- c(relmeas,1-sum(relmeas[whichO]))
          names(relmeas)[numhyp+1] <- "complement"
        }else{ #the order constrained subspaces at least partly overlap
          randomDraws <- matrix(unlist(lapply(1:numpara,function(par){
            1/rgamma(1e5,shape=shape1[par]/2,rate=scale1[par])
            #rinvgamma(1e5,shape=shape1[par]/2,scale=scale1[par])
          })),ncol=numpara)
          checksOCpost <- lapply(which(whichO),function(h){
            Rorder <- as.matrix(RrO1[[h]][,-(1+numpara)])
            if(ncol(Rorder)==1){
              Rorder <- t(Rorder)
            }
            rorder <- as.matrix(RrO1[[h]][,1+numpara])
            apply(randomDraws%*%t(Rorder) > rep(1,samsize1)%*%t(rorder),1,prod)
          })
          relmeas <- c(relmeas,sum(Reduce("+",checksOCpost) == 0) / samsize1)
          rownames(relmeas)[numhyp+1] <- "complement"
        }
      }
    }
  }
  return(relmeas)
}





