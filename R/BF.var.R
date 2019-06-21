
var_test <- function(...){ # ALLE argumenten van bartlett.test gebruiken

  temp <-  list(y = y, g = g)
  names(temp) <- c("y", "g")

  if(!is.factor(temp$g)) stop("g must be a factor")
  if(!is.numeric(temp$y)) stop("y must be a numeric")

  vars <- tapply(temp$y, temp$g, var)
  n <- table(g)
  bart <- bartlett.test(x = temp$y, g = temp$g)
  class(bart) <- c("BF_bartlett", class(bart))
  bart$vars <- vars
  bart$n <- n
  bart
}


#' @importFrom stats rchisq
#' @method BF BF_bartlett
#' @export
BF.BF_bartlett <- function(x,
                           hypothesis = NULL,
                           prior = NULL,
                           parameter = NULL,
                           ...) {

  nsim <- 1e5
  s2 <- x$vars
  n <- x$n
  J <- length(n)

  #exploratory BF for equality of variances

  #####
  #
  #FLORIAN CAN YOU ADD THE COMPUTATION OF THIS BF WHERE ALL GROUPS ARE EQUAL AGAINST UNC?
  #
  #####
  BF0u <- 3 #and the replace 3 by the actual BF

  BFtu_exploratory <- c(BF0u,1)
  names(BFtu_exploratory) <- c("homogeneity of variances","no homogeneity of variances")
  PHP_exploratory <- BFtu_exploratory / sum(BFtu_exploratory)

  if(!is.null(hypothesis)){ # execute confirmatory Bayes factor test based on hypothesis input

    parse_hyp <- parse_hypothesis(names_coef, hypothesis)
    RrList <- make_RrList2(parse_hyp)
    RrE <- RrList[[1]]
    RrO <- RrList[[2]]
    b <- 2/n
    Th <- length(RrE)
    if (is.null(prior_probs)) {
      prior_probs <- rep(1 / Th, times = Th)
    }
    logmx <- rep(NA, times = Th)
    Pr <- matrix(NA, nrow = 2, ncol = Th,
                 dimnames = list(c("Fit", "Complexity"),
                                 paste("H_", as.character(1:Th), sep = "")))
    for (h in 1:Th) {
      if (is.null(RrE[[h]])) {
        unique_vars <- as.list(1:J)
      } else {
        RrEh <- RrE[[h]][, -ncol(RrE[[h]])]
        names_coef <- colnames(RrEh)
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
      logmx[h] <- - 1 / 2 * sum((1 - unlist(blist)) * unlist(nlist)) * log(pi) +
        1 / 2 * log(prod(unlist(blist))) + sum(lgamma(df / 2) - lgamma(dfb / 2) -
                                                 1 / 2 * df * log(SS) + 1 / 2 * dfb * log(SSb))
      logmxu <- - 1 / 2 * sum((1 - b) * n) * log(pi) + 1 / 2 * log(prod(b)) +
        sum(lgamma((n - 1) / 2) - lgamma((b * n - 1) / 2) -
              1 / 2 * (n - 1) * log((n - 1) * s2) +
              1 / 2 * (b * n - 1) * log(b * (n - 1) * s2))
      if (!is.null(RrO[[h]])) {
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
        Pr[1, h] <- sum(indi_post) / nsim
        Pr[2, h] <- sum(indi_prior) / nsim
        logmx[h] <- log(Pr[1, h] / Pr[2, h]) + logmx[h]
      }
    }
    #####
    #
    # FLORIAN CAN YOU ALSO PROVIDE THE OUTPUTS 'relfit' AND 'relcomp'
    #
    ####
    #
    # FLORIAN THE COMPUTATION OF THE COMPLEMENT HYPOTHESIS CAN BE ADDED HERE.
    #
    #####
    BFtu_confirmatory <- exp(logmx - logmxu)
    names(BFtu_confirmatory) <- parse_hyp$original_hypothesis
    BFmatrix_confirmatory <- BFtu_confirmatory %*% t(1 / BFtu_confirmatory)
    diag(BFmatrix_confirmatory) <- 1
    if(is.null(prior)){
      priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
    }else{
      if(!is.numeric(prior) || length(prior)!=length(BFtu_confirmatory)){
        warning(paste0("Argument 'prior' should be numeric and of length ",as.character(length(BFtu_confirmatory)),". Equal prior probabilities are used."))
        priorprobs <- rep(1/length(BFtu_confirmatory),length(BFtu_confirmatory))
      }else{
        priorprobs <- prior
      }
    }
    PHP_confirmatory <- BFtu_confirmatory * priorprobs / sum(BFtu_confirmatory * priorprobs)
  }else{
    BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- relfit <-
      relcomp <- NULL
  }

  BFlm_out <- list(
    BFtu_exploratory=BFtu_exploratory,
    PHP_exploratory=PHP_exploratory,
    BFtu_confirmatory=BFtu_confirmatory,
    PHP_confirmatory=PHP_confirmatory,
    BFmatrix_confirmatory=BFmatrix_confirmatory,
#    relative_fit=relfit,
#    relative_complexity=relcomp,
    model=x,
    estimates=x$vars,
    hypothesis=hypothesis,
    prior=prior,
    call=match.call())

  class(BFlm_out) <- "BF"

  return(BFlm_out)
}



