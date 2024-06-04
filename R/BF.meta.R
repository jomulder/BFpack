##############################################################################
##### BF testing for meta-analysis based on Van Aert & Mulder (in prep.) #####
##### via the metafor package                                            #####
##############################################################################


#' @importFrom extraDistr rtnorm
#' @method BF rma.uni
#' @export
BF.rma.uni <- function(x,
                       hypothesis = NULL,
                       prior.hyp.explo = NULL,
                       prior.hyp.conf = NULL,
                       prior.hyp = NULL,
                       complement = TRUE,
                       log = FALSE,
                       ...){
  # x should be of class rma.uni
  if(class(x)[1]!="rma.uni"){
    stop("Only objects of class 'rma.uni' are currently supported.")
  }
  if(!is.null(x$formula.mods)){
    stop("Only an object without moderators is currently supported.")
  }
  if(!is.null(hypothesis)){
    message("Note that confirmatory testing via the 'hypothesis' argument is currently not supported for object of class 'rma.uni'.")
  }

  logIN <- log

  # check proper usage of argument 'prior.hyp.conf' and 'prior.hyp.explo'
  if(!is.null(prior.hyp.conf)){
    prior.hyp <- prior.hyp.conf
  }
  prior.hyp.explo <- process.prior.hyp.explo(prior_hyp_explo = prior.hyp.explo, model = x)

  ### Extract effect sizes and sampling variances from the metafor object
  yi <- x$yi
  vi <- x$vi

  wi <- 1/vi # Weights equal-effects model
  typ_vi <- sum(wi*(length(wi)-1))/(sum(wi)^2 - sum(wi^2)) # Typical within-study sampling variance

  ### Robbie: If statement for methods that are not "EE" or "FE"
  if (x$method != "EE" & x$method != "FE")
  {

    ### Minimal rho (add 0.001 to make sure that marginal likelihood and likelihood
    # can be evaluated)
    rho_min <- (-min(vi)+0.001)/(-min(vi)+0.001+typ_vi)

    ### Compute likelihood of model with unconstrained delta and rho = 0
    logmu0 <- marg_lik(yi = yi, vi = vi, rho = 0, rho_min = rho_min,
                       typ_vi = typ_vi)

    ### Compute likelihood of model with delta = 0 and rho unconstrained
    logm0u <- get_condpost_rho(yi = yi, vi = vi, rho_min = rho_min,
                               typ_vi = typ_vi, start_rho = x$I2/100)$logm0u

    ### Compute marginal likelihood of model with unconstrained delta and unconstrained rho
    ### Compute posterior probability of rho > 0 and rho < 0, delta unconstrained
    post_rho <- get_post_rho(yi = yi, vi = vi, rho_min = rho_min, typ_vi = typ_vi,
                             start_rho = x$I2/100)
    logmuu <- post_rho$logmuu
    prior_rho <- 1/(abs(rho_min) + 1)

    logmul <- logmuu + log(post_rho$post_rho_l/prior_rho) # Likelihood of model delta unconstrained and rho > 0
    logmus <- logmuu + log(post_rho$post_rho_s/(1-prior_rho)) # Likelihood of model delta unconstrained and rho < 0

    #get unconstrained estimates
    rhodraws <- post_rho$rhodraws
    rhostats <- c(mean(rhodraws),median(rhodraws),quantile(rhodraws,.025),quantile(rhodraws,.975))
    tau2draws <- rhodraws/(1-rhodraws)*typ_vi # Compute tau2 based on generated I^2-statistic
    mean_prior_delta <- 0
    sd_prior_delta <- sqrt(length(yi)/sum(1/(vi+mean(tau2draws))))

    mean_delta <- unlist(lapply(1:length(tau2draws), function(i){
      (mean_prior_delta/sd_prior_delta^2+sum(yi/(vi+tau2draws[i])))/
        (1/sd_prior_delta^2+sum(1/(vi+tau2draws[i])))
    }))

    sd_delta <- unlist(lapply(1:length(tau2draws), function(i){
      1/sqrt(1/sd_prior_delta^2+sum(1/(vi+tau2draws[i])))
    }))

    deltadraws <- rnorm(length(rhodraws),mean=mean_delta,sd=sd_delta)
    deltastats <- c(mean(deltadraws),median(deltadraws),quantile(deltadraws,.025),quantile(deltadraws,.975))
    uncestimates <- t(matrix(c(rhostats,deltastats),ncol=2))
    row.names(uncestimates) <- c("I^2","mu")
    colnames(uncestimates) <- c("mean","median","2.5%","97.5%")
    # posterior draws for study-specfic effects
    tau2draws <- typ_vi * rhodraws / (1 - rhodraws)
    #truncate
    tau2draws_trunc <- tau2draws[tau2draws > 0]
    deltadraws_trunc <- deltadraws[tau2draws > 0]
    numzerotau2 <- sum(tau2draws <= 0)
    deltadraws_studies <- do.call(cbind,lapply(1:length(vi),function(s){
      var_s <- 1/(1/vi[s] + 1/tau2draws_trunc)
      mean_s <- (yi[s]/vi[s] + deltadraws_trunc/tau2draws_trunc) * var_s
      c(rnorm(length(deltadraws),mean=mean_s,sd=sqrt(var_s)),rnorm(numzerotau2,mean=yi[s],sd=sqrt(vi[s])))
    }))

    ### Compute posterior probability of mu > 0 and mu < 0
    postdeltapositive <- mean(deltadraws>0)
    logmlu <- logmuu + log(postdeltapositive/.5)
    logmsu <- logmuu + log((1-postdeltapositive)/.5)

    ### Compute Bayes factors model vs. unconstrained mu and I^2
    BF0rhoUnc <- (logmu0 - logmuu)
    BF1rhoUnc <- (logmus - logmuu)
    BF2rhoUnc <- (logmul - logmuu)
    BF0deltaUnc <- (logm0u - logmuu)
    BF1deltaUnc <- (logmsu - logmuu)
    BF2deltaUnc <- (logmlu - logmuu)

    BFtu_exploratory <- matrix(c(BF0rhoUnc,BF0deltaUnc,BF1rhoUnc,BF1deltaUnc,BF2rhoUnc,BF2deltaUnc),nrow=2)
    rowmax <- apply(BFtu_exploratory,1,max)
    norm_BF_explo <- exp(BFtu_exploratory - rowmax %*% t(rep(1,ncol(BFtu_exploratory)))) *
      (rep(1,nrow(BFtu_exploratory)) %*% t(prior.hyp.explo[[1]]))

    PHP_exploratory <- norm_BF_explo / apply(norm_BF_explo,1,sum)
    colnames(BFtu_exploratory) <- c("=0","<0",">0")
    colnames(PHP_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
    row.names(BFtu_exploratory) <- row.names(PHP_exploratory) <- c("I^2","mu")
    BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- BFtable <-
      priorprobs <- hypotheses <- estimates <- NULL

    # BF_delta <- diag(rep(1, 3))
    # colnames(BF_delta) <- rownames(BF_delta) <- c("H1", "H2", "H3")
    # BF_delta[1,2] <- exp(logm0u - logmsu)
    # BF_delta[2,1] <- exp(logmsu - logm0u)
    # BF_delta[1,3] <- exp(logm0u - logmlu)
    # BF_delta[3,1] <- exp(logmlu - logm0u)
    # BF_delta[2,3] <- exp(logmsu - logmlu)
    # BF_delta[3,2] <- exp(logmlu - logmsu)

    # BF_rho <- diag(rep(1, 2))
    # colnames(BF_rho) <- rownames(BF_rho) <- c("H1", "H2")
    # BF_rho[1,2] <- exp(logmus - logmul)
    # BF_rho[2,1] <- exp(logmul - logmus)

    # sum_BF_rho <- BF1rhoUnc+BF2rhoUnc
    # PP_rho <- matrix(c(BF1rhoUnc/sum_BF_rho, BF2rhoUnc/sum_BF_rho), nrow = 1)
    # colnames(PP_rho) <- c("Pr(<0)", "Pr(>0)")

    ##############################################################################

    ### Robbie: equal-effect and fixed-effect model

  } else if (x$method == "EE" | x$method == "FE")
  {

    ### Compute likelihood of model with delta = 0 and rho = 0
    logmu <- lik(yi = yi, vi = vi, delta = 0, rho = 0, rho_min = 0, typ_vi = typ_vi)

    ### Compute likelihood of model with unconstrained delta and rho = 0
    logmu0 <- marg_lik(yi = yi, vi = vi, rho = 0, rho_min = 0,
                       typ_vi = typ_vi)

    mean_prior_delta <- 0
    sd_prior_delta <- sqrt(length(yi)/sum(1/vi))

    mean_delta <- (mean_prior_delta/sd_prior_delta^2+sum(yi/vi))/
      (1/sd_prior_delta^2+sum(1/vi))

    sd_delta <- 1/sqrt(1/sd_prior_delta^2+sum(1/vi))

    deltadraws <- rnorm(20000,mean=mean_delta,sd=sd_delta)
    deltastats <- c(mean(deltadraws),median(deltadraws),quantile(deltadraws,.025),
                    quantile(deltadraws,.975))
    uncestimates <- t(matrix(deltastats,ncol=1))
    row.names(uncestimates) <- "mu"
    colnames(uncestimates) <- c("mean","median","2.5%","97.5%")

    ### Compute posterior probability of mu > 0 and mu < 0
    postdeltapositive <- mean(deltadraws>0)
    logml <- logmu0 + log(postdeltapositive/.5)
    logms <- logmu0 + log((1-postdeltapositive)/.5)

    ### Compute Bayes factors model vs. unconstrained mu
    BF0delta <- (logmu - logmu0)
    BF1delta <- (logms - logmu0)
    BF2delta <- (logml - logmu0)

    BFtu_exploratory <- matrix(c(BF0delta,BF1delta,BF2delta),nrow=1)
    rowmax <- apply(BFtu_exploratory,1,max)
    norm_BF_explo <- exp(BFtu_exploratory - rowmax %*% t(rep(1,ncol(BFtu_exploratory)))) *
      (rep(1,nrow(BFtu_exploratory)) %*% t(prior.hyp.explo[[1]]))
    PHP_exploratory <- norm_BF_explo / apply(norm_BF_explo,1,sum)
    colnames(PHP_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
    colnames(BFtu_exploratory) <- c("H(=0)","H(<0)","H(>0)")
    row.names(BFtu_exploratory) <- row.names(PHP_exploratory) <- "mu"
    BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- BFtable <-
      priorprobs <- hypotheses <- estimates <- NULL

  }

  if(logIN==FALSE){
    BFtu_exploratory <- exp(BFtu_exploratory)
  }

  ############################

  BF_out <- list(
    BFtu_exploratory=BFtu_exploratory,
    PHP_exploratory=PHP_exploratory,
    BFtu_confirmatory=BFtu_confirmatory,
    PHP_confirmatory=PHP_confirmatory,
    BFmatrix_confirmatory=BFmatrix_confirmatory,
    BFtable_confirmatory=BFtable,
    prior.hyp.explo=prior.hyp.explo,
    prior.hyp.conf=priorprobs,
    hypotheses=hypotheses,
    estimates=uncestimates,
    model=x,
    bayesfactor="Bayes factor using uniform prior for icc & unit information prior for effect",
    parameter="between-study heterogeneity & effect size",
    log=logIN,
    call=match.call()
    #rhodraws = rhodraws,
    #deltadraws = deltadraws,
    #BF_delta = BF_delta,
    #BF_rho = BF_rho,
    #PP_rho = PP_rho
  )

  class(BF_out) <- "BF"

  BF_out

}

# rm(list = ls())

#################
### FUNCTIONS ###
#################

### Likelihood where delta is integrated out.
marg_lik <- function(yi, vi, rho, rho_min, typ_vi)
{
  tau2 <- rho/(1-rho)*typ_vi # Compute tau2 based on rho
  wi_star <- 1/(vi+tau2) # Weights random-effects model
  delta_hat <- sum(wi_star*yi)/sum(wi_star) # Estimate of delta
  k <- length(yi) # Number of studies in meta-analysis

  diagSigma <- vi+tau2
  diagSigmaInv <- 1/diagSigma

  out <- -k/2*log(2*pi) +
    -0.5*sum(log(vi+tau2)) +
    -0.5*sum((yi-delta_hat)^2 * diagSigmaInv) +
    -0.5*log(k/sum(diagSigmaInv)) +
    -log(1-rho_min) +
    -0.5*sum(diagSigmaInv)*(delta_hat^2-delta_hat^2/(1+1/k)) +
    -0.5*log(sum(diagSigmaInv)) +
    -0.5*log(1+1/k)

  return(out)
}

### Likelihood
lik <- function(yi, vi, delta, rho, rho_min, typ_vi)
{
  k <- length(yi) # Number of studies in meta-analysis
  tau2 <- rho/(1-rho)*typ_vi # Compute tau2 based on rho

  out <- -(k+1)/2*log(2*pi) +
    -0.5*sum(log(prod(vi+tau2))) +
    -0.5*sum((yi-delta)^2/(vi+tau2)) +
    -log(1-rho_min) +
    -0.5*log(k/sum((vi+tau2)^(-1))) +
    -0.5*delta^2/(k/sum((vi+tau2)^(-1)))

  return(out)
}

### General method-of-moments estimate (Eq. 6 in DerSimonian and Kacker, 2007)
MM <- function(yi, vi, ai)
{
  yw <- sum(ai*yi)/sum(ai)
  tau2_hat <- (sum(ai*(yi-yw)^2)-(sum(ai*vi)-sum(ai^2*vi)/sum(ai)))/(sum(ai)-sum(ai^2)/sum(ai))
  return(tau2_hat)
}


get_post_rho <- function(yi, vi, rho_min, typ_vi, start_rho, iters = 20000)
{

  rho_s <- numeric(iters) # Empty object for storing results
  check1 <- 100
  #burn in
  itersBI <- 5000 # length of burn-in
  rhoBI <- start_rho # some starting value for rho...
  sdstep <- .1
  acceptMat <- rep(0, length = iters + itersBI)
  rhoBIseq <- rep(0, length = itersBI)
  sdstepseq <- rep(0, length = (iters + itersBI) / check1)
  sdsteptel <- 1
  upper1 <- .5 # define region where the sdstep does not have to be changed
  lower1 <- .15
  #start burn-in
  for(i in 1:itersBI){
    #draw candidate from truncated normal
    rho_star <- rtnorm(1, mean = rhoBI, sd = sdstep, a = rho_min, b = 1)
    #evaluate Metropolis-Hastings acceptance probability
    R_MH <- exp( marg_lik(yi = yi, vi = vi, rho = rho_star,
                          rho_min = rho_min, typ_vi = typ_vi) -
                   marg_lik(yi = yi, vi = vi, rho = rhoBI,
                            rho_min = rho_min, typ_vi = typ_vi) ) *
      (pnorm(1,mean=rhoBI,sd=sdstep) - pnorm(rho_min,mean=rhoBI,sd=sdstep)) /
      (pnorm(1,mean=rho_star,sd=sdstep) - pnorm(rho_min,mean=rho_star,sd=sdstep))
    rhoBI <- ifelse(runif(1) < R_MH, rho_star, rhoBI)
    acceptMat[i] <- rho_star == rhoBI
    rhoBIseq[i] <- rhoBI
    #if needed update random walk sd depending on acceptance rate
    if(ceiling(i/check1)==i/check1){
      probs <- mean(acceptMat[(i-check1+1):i])
      if(probs>upper1){
        sdstep <- sdstep * ( (probs-upper1)/(1-upper1) + 1)
      }else if(probs < lower1){
        sdstep <- sdstep * 1 / ( 2 - probs/lower1 )
      }
      sdstep <- ifelse(sdstep>1,1,sdstep)
      sdstepseq[sdsteptel] <- sdstep
      sdsteptel <- sdsteptel + 1
    }
  }
  #now actual drawing
  rho_s[1] <- rhoBI
  for(i in 2:iters){
    #draw candidate from truncated normal
    rho_star <- rtnorm(1, mean = rho_s[i-1], sd = sdstep, a = rho_min, b = 1)
    #evaluate Metropolis-Hastings acceptance probability
    R_MH <- exp( marg_lik(yi = yi, vi = vi, rho = rho_star,
                          rho_min = rho_min, typ_vi = typ_vi) -
                   marg_lik(yi = yi, vi = vi, rho = rho_s[i-1],
                            rho_min = rho_min, typ_vi = typ_vi) ) *
      (pnorm(1,mean=rho_s[i-1],sd=sdstep) - pnorm(rho_min,mean=rho_s[i-1],sd=sdstep)) /
      (pnorm(1,mean=rho_star,sd=sdstep) - pnorm(rho_min,mean=rho_star,sd=sdstep))
    rho_s[i] <- ifelse(runif(1)<R_MH,rho_star,rho_s[i-1])
    acceptMat[i] <- rho_star==rho_s[i]
    #if needed update random walk sd depending on acceptance rate
    if(ceiling(i/check1)==i/check1){
      probs <- mean(acceptMat[(i-check1+1):i])
      if(probs>upper1){
        sdstep <- sdstep * ( (probs-upper1)/(1-upper1) + 1)
      }else if(probs < lower1){
        sdstep <- sdstep * 1/( 2 - probs/lower1 )
      }
      #given the bounds on rho ststep should not have to be larger than 1
      #this also avoids unlimited growth of sdstep in case of relatively diffuse distribution
      sdstep <- ifelse(sdstep > 1,1,sdstep)
      sdstepseq[sdsteptel] <- sdstep
      sdsteptel <- sdsteptel + 1
    }
  }
  ### Compute posterior model probabilities
  post_rho_s <- mean(rho_s < 0)
  post_rho_l <- 1 - post_rho_s

  ### Compute marginal likelihood of rho and delta unconstrained
  meanICC <- mean(rho_s)
  varICC <- var(rho_s)
  shape1IM <- - ((rho_min-meanICC)*((meanICC-1)*meanICC-meanICC*rho_min+rho_min+varICC)/
                   ((rho_min-1)*varICC))
  shape2IM <- - shape1IM*(meanICC-1)/(meanICC-rho_min)
  factor1 <- .6
  shape1IM <- shape1IM * factor1
  shape2IM <- shape2IM * factor1
  ISnum <- 1e4
  #importance sampler from stretched beta distribution
  ISdraws <- (rbeta(ISnum,shape1IM,shape2IM) * (1 - rho_min) + rho_min) * .99999
  logintegrands <- unlist(lapply(1:ISnum,function(s){

    marg_lik(yi = yi, vi = vi, rho = ISdraws[s], rho_min = rho_min, typ_vi = typ_vi) -
      dbeta((ISdraws[s]-rho_min)/(1-rho_min),shape1=shape1IM,shape2=shape2IM,log=TRUE) +
      log(1/(1-rho_min))
  }))
  #Sys.time() - wer
  logmuu <- log(mean(exp(logintegrands-max(logintegrands)))) + max(logintegrands)

  return(list(post_rho_s = post_rho_s, post_rho_l = post_rho_l, rhodraws = rho_s, logmuu = logmuu))
}

get_condpost_rho <- function(yi, vi, rho_min, typ_vi, start_rho, iters = 20000)
{
  #condition on delta = 0
  rho_s <- numeric(iters) # Empty object for storing results
  check1 <- 100
  #burn in
  itersBI <- 5000 # length of burn-in
  rhoBI <- start_rho # some starting value for rho...
  sdstep <- .1
  acceptMat <- rep(0, length = iters + itersBI)
  rhoBIseq <- rep(0, length = itersBI)
  sdstepseq <- rep(0, length = (iters + itersBI) / check1)
  sdsteptel <- 1
  upper1 <- .5 # define region where the sdstep does not have to be changed
  lower1 <- .15
  #start burn-in
  for(i in 1:itersBI){
    #draw candidate from truncated normal
    rho_star <- rtnorm(1, mean = rhoBI, sd = sdstep, a = rho_min, b = 1)
    #evaluate Metropolis-Hastings acceptance probability
    R_MH <- exp(lik(yi = yi, vi = vi, delta = 0, rho = rho_star,
                    rho_min = rho_min, typ_vi = typ_vi) -
                  lik(yi = yi, vi = vi, delta = 0, rho = rhoBI,
                      rho_min = rho_min, typ_vi = typ_vi) ) *
      (pnorm(1,mean=rhoBI,sd=sdstep) - pnorm(rho_min,mean=rhoBI,sd=sdstep)) /
      (pnorm(1,mean=rho_star,sd=sdstep) - pnorm(rho_min,mean=rho_star,sd=sdstep))
    rhoBI <- ifelse(runif(1) < R_MH, rho_star, rhoBI)
    acceptMat[i] <- rho_star == rhoBI
    rhoBIseq[i] <- rhoBI
    #if needed update random walk sd depending on acceptance rate
    if(ceiling(i/check1)==i/check1){
      probs <- mean(acceptMat[(i-check1+1):i])
      if(probs>upper1){
        sdstep <- sdstep * ( (probs-upper1)/(1-upper1) + 1)
      }else if(probs < lower1){
        sdstep <- sdstep * 1 / ( 2 - probs/lower1 )
      }
      sdstep <- ifelse(sdstep>1,1,sdstep)
      sdstepseq[sdsteptel] <- sdstep
      sdsteptel <- sdsteptel + 1
    }
  }
  #now actual drawing
  rho_s[1] <- rhoBI
  for(i in 2:iters){
    #draw candidate from truncated normal
    rho_star <- rtnorm(1, mean = rho_s[i-1], sd = sdstep, a = rho_min, b = 1)
    #evaluate Metropolis-Hastings acceptance probability
    R_MH <- exp( lik(yi = yi, vi = vi, delta = 0, rho = rho_star,
                     rho_min = rho_min, typ_vi = typ_vi) -
                   lik(yi = yi, vi = vi, delta = 0, rho = rho_s[i-1],
                       rho_min = rho_min, typ_vi = typ_vi) ) *
      (pnorm(1,mean=rho_s[i-1],sd=sdstep) - pnorm(rho_min,mean=rho_s[i-1],sd=sdstep)) /
      (pnorm(1,mean=rho_star,sd=sdstep) - pnorm(rho_min,mean=rho_star,sd=sdstep))
    rho_s[i] <- ifelse(runif(1) < R_MH, rho_star, rho_s[i-1])
    acceptMat[i] <- rho_star==rho_s[i]
    #if needed update random walk sd depending on acceptance rate
    if(ceiling(i/check1)==i/check1){
      probs <- mean(acceptMat[(i-check1+1):i])
      if(probs>upper1){
        sdstep <- sdstep * ( (probs-upper1)/(1-upper1) + 1)
      }else if(probs < lower1){
        sdstep <- sdstep * 1/( 2 - probs/lower1 )
      }
      #given the bounds on rho ststep should not have to be larger than 1
      #this also avoids unlimited growth of sdstep in case of relatively diffuse distribution
      sdstep <- ifelse(sdstep > 1,1,sdstep)
      sdstepseq[sdsteptel] <- sdstep
      sdsteptel <- sdsteptel + 1
    }
  }

  ### Compute marginal likelihood of rho and delta unconstrained
  meanICC <- mean(rho_s)
  varICC <- var(rho_s)
  shape1IM <- - ((rho_min-meanICC)*((meanICC-1)*meanICC-meanICC*rho_min+rho_min+varICC)/
                   ((rho_min-1)*varICC))
  shape2IM <- - shape1IM*(meanICC-1)/(meanICC-rho_min)
  factor1 <- .6
  shape1IM <- shape1IM * factor1
  shape2IM <- shape2IM * factor1
  ISnum <- 1e4
  #importance sampler from stretched beta distribution
  ISdraws <- (rbeta(ISnum,shape1IM,shape2IM) * (1 - rho_min) + rho_min) * .99999
  logintegrands <- unlist(lapply(1:ISnum,function(s){
    lik(yi = yi, vi = vi, delta = 0, rho = ISdraws[s], rho_min = rho_min, typ_vi = typ_vi) -
      dbeta((ISdraws[s]-rho_min)/(1-rho_min),shape1=shape1IM,shape2=shape2IM,log=TRUE) +
      log(1/(1-rho_min))
  }))
  #Sys.time() - wer
  logm0u <- log(mean(exp(logintegrands-max(logintegrands)))) + max(logintegrands)

  return(list(rhodraws = rho_s, logm0u = logm0u))
}
