
#' @importFrom metaBMA prior
#' @importFrom stats dlogis rlnorm dlnorm
#' @importFrom extraDistr rtnorm
#' @importFrom MASS fitdistr
#' @importFrom methods is
#' @describeIn BF BF S3 method for an object of class 'rma.uni'
#' @method BF rma.uni
#' @export
BF.rma.uni <- function(x,
                    hypothesis = NULL,
                    prior.hyp.explo = NULL,
                    prior.hyp.conf = NULL,
                    prior.hyp = NULL,
                    complement = TRUE,
                    log = FALSE,
                    BF.type,
                    iter = 2e4,
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

  if(!exists("BF.type")){
    stop("The argument 'BF.type' is missing. See documentation. See ?BF")
  }

  logIN <- log

  # check proper usage of argument 'prior.hyp.conf' and 'prior.hyp.explo'
  if(!is.null(prior.hyp.conf)){
    prior.hyp <- prior.hyp.conf
  }
  prior.hyp.explo <- process.prior.hyp.explo(prior_hyp_explo = prior.hyp.explo, model = x)

  ### Extract effect sizes and sampling variances from the metafor object
  yi <- c(x$yi)
  vi <- x$vi
  ni <- x$ni #only used in case of unit-information prior: BF.type = "unit.info"
  N <- sum(ni)
  K <- length(vi)

  wi <- 1/vi # Weights equal-effects model
  typ_vi <- sum(wi*(length(wi)-1))/(sum(wi)^2 - sum(wi^2)) # Typical within-study sampling variance
  est.mu <- sum(yi*wi) / sum(wi)
  se.mu <- x$se

  est.tau2 <- x$tau2
  se.tau2 <- x$se.tau2 # check if we can use it as step sd for posterior sampling of tau
  if(is.na(se.tau2)){se.tau2 <- .5}
  tau2.min <- -min(vi)

  # prior tau2
  prior.tau2.Jeffreys1 <- function(tau2.arg,log=FALSE){
    ptau2 <- -1/K*sum(log(vi+tau2.arg))
    #ptau2 <- .5*log(sum(1/(vi+tau2.arg)^2)) # Tibshirani (1989)
    #ptau2 <- 0 #uniform prior
    ifelse(log,ptau2,exp(ptau2))
    } # Berger & Deely (1988)


  if(is(BF.type,"character")){
    if(BF.type == "stand.effect"){
      # conjugate normal prior assuming average effects of about .5
      prior.mu <- function(x,tau2,log=FALSE){dnorm(x, mean = 0, sd = .5, log = log)}
      prior.muGR0 <- .5
      bayesfactor.name <- "Bayes factor based on a normal prior (mu ~ norm(mean=0, sd=0.5))"
      parameter.name <- "standardized effect"
    }else if(BF.type == "log.odds"){
      # Student t prior which approximates the implied distribution of the log odds ratio based on uniform success probabilities
      prior.mu <- function(x,tau2,log=FALSE){
        dt1 <- dt(x/2.36,df=13.1,log=log)-log(2.36)
        ifelse(log,dt1,exp(dt1))
        }
      prior.muGR0 <- .5
      bayesfactor.name <- "Bayes factor based on uniform priors for success probabilities (log.odds ~ t(0,2.36,13.1))."
      parameter.name <- "log odds"
    }else if(BF.type == "correlation"){
      # logistic prior for the Fisher transformed correlation corresponding to a uniform prior for the correlation in (-1,1)
      prior.mu <- function(x,tau2,log=FALSE){dlogis(x, scale = .5, log = log)}
      prior.muGR0 <- .5
      bayesfactor.name <- "Bayes factor based on a uniform prior for the correlation in (-1,1) (Fisher(cor)~logis(0.5))."
      parameter.name <- "correlation"
    }else if(BF.type == "unit.info"){
      if(N==0){
        stop("In order to use a unit-information prior, the sample sizes need to be part of the element 'ni'
              of the 'rma.uni' object (e.g., using the 'ni' argument of the 'rma' function).")
      }
      prior.mu <- function(x,tau2,log=FALSE){dnorm(x, mean = 0, sd = sqrt(N/sum(1/(vi+tau2))), log = log)}
      prior.muGR0 <- .5
      bayesfactor.name <- "Bayes factor based on unit information prior."
      parameter.name <- "general parameter"
    }else{
      stop("The argument 'BF.type' is not correctly specified for an object of type 'rma.uni'. See documentation. ?BF")
    }
  }else if(is(BF.type,"prior")){
    prior.mu <- function(x,tau2,log=FALSE){BF.type(x,log)}
    seq1 <- seq(0,1e3,length=1e5)
    seq1 <- seq1+(seq1[2]-seq1[1])/2
    prior.muGR0 <- sum(prior.mu(seq1,0))*(seq1[2]-seq1[1]) #riemann estimate of prior prob that mu > 0.
    bayesfactor.name <- paste0("Bayes factor based on manually chosen prior")
    parameter.name <- "general parameter"
    message("Be sure that the specified prior is not truncated in a specific interval.")
  }else{
    stop("The argument 'BF.type' is not correctly specified for an object of type 'rma.uni'. See documentation. ?BF")
  }

  #exploratory testing
  if (x$method != "EE" & x$method != "FE") ### random effects / marema model
  {

    #### MAREMA MODEL

    # initial values for MCMC run
    start_mu <- est.mu
    sdstep.mu <- x$se
    start_tau2 <- ifelse(x$tau2>0, x$tau2,.05)
    sdstep.tau2 <- x$se.tau2

    # marginal likelihood full unconstrained model. marema.unc.post.draws
    marg.like.marema.full <- log_marg_like_full(yi, vi, tau2_min=tau2.min,
                                                prior.mu=prior.mu,
                                                prior.tau2=prior.tau2.Jeffreys1,
                                                start_mu=start_mu, start_tau2=start_tau2,
                                                sdstep.mu=sdstep.mu, sdstep.tau2=sdstep.tau2,
                                                burnin = round(iter/2), iters1=iter, iters2=iter)

    # marginal likelihood under mu = 0
    marg.like.marema.muEQ0 <- log_marg_like_cond.mu(yi, vi, tau2_min=tau2.min,
                                                    muIN=0, prior.tau2=prior.tau2.Jeffreys1,
                                                    start_tau2=start_tau2+start_mu^2, sdstep.tau2=sdstep.tau2,
                                                    burnin=iter, iters1=iter, iters2=iter)

    # marginal likelihoods when constraining mu
    # mu = 0
    marg.like.marema.muEQ0 <- marg.like.marema.muEQ0[[1]]
    # mu > 0
    marg.like.marema.muGR0 <- marg.like.marema.full[[1]] + log(mean(marg.like.marema.full$post.draws[,1]>0)/prior.muGR0)
    # mu < 0
    marg.like.marema.muSM0 <- marg.like.marema.full[[1]] + log(mean(marg.like.marema.full$post.draws[,1]<0)/(1-prior.muGR0))

    #### random effects MODEL

    # marginal likelihood full unconstrained model
    marg.like.ranef.full <- log_marg_like_full(yi, vi, tau2_min=0,
                                               prior.mu = prior.mu,
                                               prior.tau2 = prior.tau2.Jeffreys1,
                                               start_mu = start_mu, start_tau2 = ifelse(start_tau2<0,sdstep.tau2,start_tau2),
                                               sdstep.mu = sdstep.mu, sdstep.tau2 = sdstep.tau2,
                                               burnin=round(iter/2), iters1=iter, iters2=iter)

    # marginal likelihood under mu = 0
    marg.like.ranef.muEQ0 <- log_marg_like_cond.mu(yi, vi, tau2_min=0,
                                                   muIN=0, prior.tau2=prior.tau2.Jeffreys1,
                                                   start_tau2=ifelse(start_tau2<0,sdstep.tau2,start_tau2)+start_mu^2, sdstep.tau2=sdstep.tau2,
                                                   burnin = round(iter/2), iters1=iter, iters2=iter)

    marg.like.ranef.muEQ0 <- marg.like.ranef.muEQ0[[1]]
    # mu > 0
    marg.like.ranef.muGR0 <- marg.like.ranef.full[[1]] + log(mean(marg.like.ranef.full$post.draws[,1]>0)/prior.muGR0)
    # mu < 0
    marg.like.ranef.muSM0 <- marg.like.ranef.full[[1]] + log(mean(marg.like.ranef.full$post.draws[,1]<0)/(1-prior.muGR0))

    ### Compute Bayes factors model vs. unconstrained mu and I^2 (marema)
    BF0u_ranef <- marg.like.ranef.muEQ0 - marg.like.ranef.full[[1]]
    BF1u_ranef <- marg.like.ranef.muSM0 - marg.like.ranef.full[[1]]
    BF2u_ranef <- marg.like.ranef.muGR0 - marg.like.ranef.full[[1]]
    BF0u_marema <- marg.like.marema.muEQ0 - marg.like.marema.full[[1]]
    BF1u_marema <- marg.like.marema.muSM0 - marg.like.marema.full[[1]]
    BF2u_marema <- marg.like.marema.muGR0 - marg.like.marema.full[[1]]

    BFtu_exploratory_marema <- matrix(c(BF0u_marema,BF1u_marema,BF2u_marema),nrow=1)
    rowmax <- apply(BFtu_exploratory_marema,1,max)
    norm_BF_explo <- exp(BFtu_exploratory_marema - rowmax %*% t(rep(1,ncol(BFtu_exploratory_marema)))) *
      (rep(1,nrow(BFtu_exploratory_marema)) %*% t(prior.hyp.explo[[1]]))
    PHP_exploratory_marema <- norm_BF_explo / apply(norm_BF_explo,1,sum)

    BFtu_exploratory_ranef <- matrix(c(BF0u_ranef,BF1u_ranef,BF2u_ranef),nrow=1)
    rowmax <- apply(BFtu_exploratory_ranef,1,max)
    norm_BF_explo <- exp(BFtu_exploratory_ranef - rowmax %*% t(rep(1,ncol(BFtu_exploratory_ranef)))) *
      (rep(1,nrow(BFtu_exploratory_ranef)) %*% t(prior.hyp.explo[[1]]))
    PHP_exploratory_ranef <- norm_BF_explo / apply(norm_BF_explo,1,sum)

    BFtu_exploratory <- rbind(BFtu_exploratory_marema,BFtu_exploratory_ranef)
    PHP_exploratory <- rbind(PHP_exploratory_marema,PHP_exploratory_ranef)

    colnames(BFtu_exploratory) <- c("=0","<0",">0")
    colnames(PHP_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
    row.names(BFtu_exploratory) <- row.names(PHP_exploratory) <- c("mu (marema)", "mu (ranef)")

    if(logIN == FALSE){
      BFtu_exploratory <- exp(BFtu_exploratory)
    }

    #store descriptives of estimates
    estimates_marema <- cbind(apply(marg.like.marema.full$post.draws,2,mean),
                              apply(marg.like.marema.full$post.draws,2,median),
                              apply(marg.like.marema.full$post.draws,2,quantile,.025),
                              apply(marg.like.marema.full$post.draws,2,quantile,.975),
                              apply(marg.like.marema.full$post.draws,2,function(x){mean(x>0)}))
    row.names(estimates_marema) <- c("mu      (marema)", "tau2    (marema)")
    estimates_ranef <- cbind(apply(marg.like.ranef.full$post.draws,2,mean),
                             apply(marg.like.ranef.full$post.draws,2,median),
                             apply(marg.like.ranef.full$post.draws,2,quantile,.025),
                             apply(marg.like.ranef.full$post.draws,2,quantile,.975),
                             apply(marg.like.ranef.full$post.draws,2,function(x){mean(x>0)}))
    row.names(estimates_ranef) <- c("mu      (ranef)", "tau2    (ranef)")
    colnames(estimates_marema) <- colnames(estimates_ranef) <- c("mean","median","2.5%","97.5%","Pr(>0)")
    uncestimates <- estimates_marema
    # estimates of separate group means
    #under marema model
    marema_post.draws_mu <- marg.like.marema.full$post.draws[,1]
    marema_post.draws_tau2 <- marg.like.marema.full$post.draws[,2]
    tau2draws_trunc <- marg.like.marema.full$post.draws[marema_post.draws_tau2>0,2]
    mudraws_trunc <- marg.like.marema.full$post.draws[marema_post.draws_tau2>0,1]
    numzerotau2 <- sum(marema_post.draws_tau2 <= 0)
    mudraws_studies <- do.call(cbind,lapply(1:length(vi),function(s){
      var_s <- 1/(1/vi[s] + 1/tau2draws_trunc)
      mean_s <- (yi[s]/vi[s] + mudraws_trunc/tau2draws_trunc) * var_s
      draws_s <- c(rnorm(length(mudraws_trunc),mean=mean_s,sd=sqrt(var_s)),
                   marema_post.draws_mu[marema_post.draws_tau2<0])
      c(draws_s,mean(draws_s),median(draws_s),quantile(draws_s,.025),quantile(draws_s,.975),mean(draws_s>0))
    }))
    uncestimates <- rbind(uncestimates,t(mudraws_studies[length(marema_post.draws_mu)+1:5,]))
    row.names(uncestimates)[-(1:2)] <- paste0("theta_",1:ncol(mudraws_studies)," (marema)")
    #under random effects model
    uncestimates <- rbind(uncestimates,estimates_ranef)
    ranef_post.draws_mu <- marg.like.ranef.full$post.draws[,1]
    ranef_post.draws_tau2 <- marg.like.ranef.full$post.draws[,2]
    mudraws_studies <- do.call(cbind,lapply(1:length(vi),function(s){
      var_s <- 1/(1/vi[s] + 1/ranef_post.draws_tau2)
      mean_s <- (yi[s]/vi[s] + ranef_post.draws_mu/ranef_post.draws_tau2) * var_s
      draws_s <- rnorm(length(ranef_post.draws_mu),mean=mean_s,sd=sqrt(var_s))
      c(draws_s,mean(draws_s),median(draws_s),quantile(draws_s,.025),quantile(draws_s,.975),mean(draws_s>0))
    }))
    uncestimates <- rbind(uncestimates,t(mudraws_studies[length(ranef_post.draws_mu)+1:5,]))
    row.names(uncestimates)[-(1:(4+length(vi)))] <- paste0("theta_",1:ncol(mudraws_studies)," (ranef)")

    #store posterior draws based on noninformative improper prior
    post.draws <- list()
    post.draws$marema <- marg.like.marema.full$post.draws
    post.draws$ranef <- marg.like.ranef.full$post.draws

  } else if (x$method == "EE" | x$method == "FE") { ### common effect model

    # testing
    marg.like.comef.muEQ0 <- log_marg_like_cond.tau2(yi, vi, tau2IN=0,
                                                     prior.mu,
                                                     start_mu=est.mu, sdstep.mu=se.mu,
                                                     burnin = round(iter/2), iters1 = iter, iters2 = iter)

    ### log marginal likelihoods
    # mu = 0
    marg.like.muEQ0 <- sum(dnorm(yi,mean=0,sd=sqrt(vi+0),log=TRUE))
    # mu > 0
    marg.like.muGR0 <- marg.like.comef.muEQ0[[1]] + log(mean(marg.like.comef.muEQ0[[2]]>0)/prior.muGR0)
    # mu < 0
    marg.like.muSM0 <- marg.like.comef.muEQ0[[1]] + log(mean(marg.like.comef.muEQ0[[2]]<0)/(1-prior.muGR0))

    ### Compute Bayes factors model vs. unconstrained model
    BF0u_comef <- marg.like.muEQ0 - marg.like.comef.muEQ0[[1]]
    BF1u_comef <- marg.like.muSM0 - marg.like.comef.muEQ0[[1]]
    BF2u_comef <- marg.like.muGR0 - marg.like.comef.muEQ0[[1]]

    BFtu_exploratory <- matrix(c(BF0u_comef,BF1u_comef,BF2u_comef),nrow=1)
    rowmax <- apply(BFtu_exploratory,1,max)
    norm_BF_explo <- exp(BFtu_exploratory - rowmax %*% t(rep(1,ncol(BFtu_exploratory)))) *
      (rep(1,nrow(BFtu_exploratory)) %*% t(prior.hyp.explo[[1]]))
    PHP_exploratory <- norm_BF_explo / apply(norm_BF_explo,1,sum)
    colnames(BFtu_exploratory) <- c("=0","<0",">0")
    colnames(PHP_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
    row.names(BFtu_exploratory) <- row.names(PHP_exploratory) <- c("mu (comef)")

    #descriptives of estimates
    uncestimates <- cbind(mean(marg.like.comef.muEQ0$post.draws), median(marg.like.comef.muEQ0$post.draws),
                          quantile(marg.like.comef.muEQ0$post.draws,.025),quantile(marg.like.comef.muEQ0$post.draws,.975),
                          mean(marg.like.comef.muEQ0$post.draws>0))
    colnames(uncestimates) <- c("mean","median","2.5%","97.5%","Pr(>0)")
    row.names(uncestimates) <- "mu"

    #store posterior draws based on noninformative improper prior
    post.draws <- marg.like.comef.muEQ0$post.draws
    colnames(post.draws) <- row.names(uncestimates)

    if(logIN == FALSE){
      BFtu_exploratory <- exp(BFtu_exploratory)
    }

  }

  if(!is.null(hypothesis)){
    message("The 'hypothesis' argument for this model class is not yet supported. Please inform us about the types of hypotheses you wish to test :)")
  }

  BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- BFtable <-
    priorprobs <- hypotheses <- NULL

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
    bayesfactor=bayesfactor.name,
    parameter=parameter.name,
    log=logIN,
    call=match.call(),
    post.draws=post.draws
  )

  class(BF_out) <- "BF"

  BF_out

}

#1. unconstrained sampling under marema using 'prior.mu' as prior for mu and a noninformative prior for tau in (prior.tau2,Inf)
gibbs_unc_prior <- function(yi, vi, tau2_min,
                            prior.mu, prior.tau2,
                            start_mu, start_tau2,
                            sdstep.mu, sdstep.tau2,
                            burnin, iters){

  #check sd for random walk every 100 iterations
  check1 <- 100
  upper1 <- .5 # define region where the sdstep does not have to be changed
  lower1 <- .15
  sdsteptel <- 1
  #stretched beta prior for rho is a uniform: alpha1 = 1 & beta1 = 1

  #initialization
  mu <- start_mu
  if(start_tau2>tau2_min){
    tau2 <- start_tau2
  }else{
    tau2 <- sdstep.tau2
  }
  cur_vars <- vi+tau2
  acceptMat <- matrix(0, nrow = iters + burnin, ncol = 2)
  sdstepseq <- matrix(0, nrow = (iters + burnin) / check1, 2)
  sdstep <- c(sdstep.mu, sdstep.tau2)
  store <- matrix(NA, ncol=2, nrow = iters + burnin)
  colnames(store) <- colnames(acceptMat) <- colnames(sdstepseq) <- c("mu","tau2")

  for(i in 1:(burnin+iters)){
    #draw candidate tau2 from truncated normal
    tau2_star <- rtnorm(1, mean = tau2, sd = sdstep[2], a = tau2_min)

    #evaluate Metropolis-Hastings acceptance probability
    can_vars <- vi+tau2_star
    R_MH <- exp( sum(dnorm(yi,mean=mu,sd=sqrt(can_vars),log=TRUE)) -
                   sum(dnorm(yi,mean=mu,sd=sqrt(cur_vars),log=TRUE)) +
                   prior.mu(mu,tau2_star,log=TRUE) - prior.mu(mu,tau2,log=TRUE) +
                   prior.tau2(tau2_star,log=TRUE) - prior.tau2(tau2,log=TRUE)) *
      (1 - pnorm(tau2_min,mean=tau2,sd=sdstep[2])) / (1 - pnorm(tau2_min,mean=tau2_star,sd=sdstep[2]))
    tau2 <- ifelse(runif(1) < R_MH, tau2_star, tau2)
    acceptMat[i,2] <- tau2_star == tau2
    cur_vars <- vi+tau2

    #draw mu
    mu_star <- rnorm(1, mean = mu, sd = sdstep[1])
    R_MH <- exp( sum(dnorm(yi,mean=mu_star,sd=sqrt(cur_vars),log=TRUE)) -
                   sum(dnorm(yi,mean=mu,sd=sqrt(cur_vars),log=TRUE)) +
      prior.mu(mu_star,tau2,log=TRUE) - prior.mu(mu,tau2,log=TRUE)
      )
    mu <- ifelse(runif(1) < R_MH, mu_star, mu)
    acceptMat[i,1] <- mu_star == mu

    #if needed update random walk sd depending on acceptance rate
    if(ceiling(i/check1)==i/check1){
      probs <- apply(as.matrix(acceptMat[(i-check1+1):i,]),2,mean)
      sdstep[probs>upper1] <- sdstep[probs>upper1] * ( (probs[probs>upper1]-upper1)/(1-upper1) + 1)
      sdstep[probs<lower1] <- sdstep[probs<lower1] * 1/( 2 - (probs[probs<lower1])/lower1 )
      # store sdstep
      sdstepseq[sdsteptel,] <- sdstep
      sdsteptel <- sdsteptel + 1
    }

    store[i,] <- c(mu,tau2)

  }

  return(list(postdraws = as.matrix(store[-(1:burnin),]), sdstepseq=sdstepseq))

}

#2. unconstrained sampling under marema while fixing 'mu' and a noninformative prior for tau in (prior.tau2,Inf)
gibbs_cond.mu_prior <- function(yi, vi, tau2_min,
                                muIN, prior.tau2,
                                start_tau2, sdstep.tau2, burnin, iters){

  #check sd for random walk every 100 iterations
  check1 <- 100
  upper1 <- .5 # define region where the sdstep does not have to be changed
  lower1 <- .15
  sdsteptel <- 1
  #stretched beta prior for rho is a uniform: alpha1 = 1 & beta1 = 1

  #initialization
  mu <- muIN
  if(start_tau2>tau2_min){
    tau2 <- start_tau2
  }else{
    tau2 <- sdstep.tau2
  }
  cur_vars <- vi+tau2
  acceptMat <- matrix(0, nrow = iters + burnin, ncol = 1)
  sdstepseq <- matrix(0, nrow = (iters + burnin) / check1, 1)
  sdstep <- sdstep.tau2
  store <- matrix(NA, ncol=1, nrow = iters + burnin)
  colnames(store) <- colnames(acceptMat) <- colnames(sdstepseq) <- c("tau2")

  for(i in 1:(burnin+iters)){
    #draw rho
    #draw candidate from truncated normal
    tau2_star <- rtnorm(1, mean = tau2, sd = sdstep, a = tau2_min)

    #evaluate Metropolis-Hastings acceptance probability
    can_vars <- vi+tau2_star
    R_MH <- exp( sum(dnorm(yi,mean=mu,sd=sqrt(can_vars),log=TRUE)) -
                   sum(dnorm(yi,mean=mu,sd=sqrt(cur_vars),log=TRUE)) ) *
      (1 - pnorm(tau2_min,mean=tau2,sd=sdstep)) / (1 - pnorm(tau2_min,mean=tau2_star,sd=sdstep)) *
      prior.tau2(tau2_star) / prior.tau2(tau2)
    tau2 <- ifelse(runif(1) < R_MH, tau2_star, tau2)
    acceptMat[i,1] <- tau2_star == tau2
    cur_vars <- vi+tau2

    #if needed update random walk sd depending on acceptance rate
    if(ceiling(i/check1)==i/check1){
      probs <- mean(as.matrix(acceptMat[(i-check1+1):i,]))
      sdstep[probs>upper1] <- sdstep[probs>upper1] * ( (probs[probs>upper1]-upper1)/(1-upper1) + 1)
      sdstep[probs<lower1] <- sdstep[probs<lower1] * 1/( 2 - (probs[probs<lower1])/lower1 )
      # store sdstep
      sdstepseq[sdsteptel,1] <- sdstep
      sdsteptel <- sdsteptel + 1
    }

    store[i,] <- c(tau2)

  }

  return(list(postdraws = as.matrix(store[-(1:burnin),]), sdstepseq=sdstepseq))

}

#3. unconstrained sampling under marema using 'prior.mu' as prior while fixing 'tau2'
gibbs_cond.tau_prior <- function(yi, vi, tau2IN,
                                 prior.mu,
                                 start_mu, sdstep.mu, burnin, iters){

  #check sd for random walk every 100 iterations
  check1 <- 100
  upper1 <- .5 # define region where the sdstep does not have to be changed
  lower1 <- .15
  sdsteptel <- 1
  #stretched beta prior for rho is a uniform: alpha1 = 1 & beta1 = 1

  #initialization
  mu <- start_mu
  tau2 <- tau2IN
  cur_vars <- vi+tau2
  acceptMat <- matrix(0, nrow = iters + burnin, ncol = 1)
  sdstepseq <- matrix(0, nrow = (iters + burnin) / check1, 1)
  sdstep <- c(sdstep.mu)
  store <- matrix(NA, ncol=1, nrow = iters + burnin)
  colnames(store) <- colnames(acceptMat) <- colnames(sdstepseq) <- c("mu")

  for(i in 1:(burnin+iters)){

    #draw mu
    mu_star <- rnorm(1, mean = mu, sd = sdstep[1])

    R_MH <- exp( sum(dnorm(yi,mean=mu_star,sd=sqrt(cur_vars),log=TRUE)) -
                   sum(dnorm(yi,mean=mu,sd=sqrt(cur_vars),log=TRUE)) ) *
      prior.mu(mu_star,tau2) / prior.mu(mu,tau2)
    mu <- ifelse(runif(1) < R_MH, mu_star, mu)
    acceptMat[i,1] <- mu_star == mu

    #if needed update random walk sd depending on acceptance rate
    if(ceiling(i/check1)==i/check1){
      probs <- mean(as.matrix(acceptMat[(i-check1+1):i,]))
      sdstep[probs>upper1] <- sdstep[probs>upper1] * ( (probs[probs>upper1]-upper1)/(1-upper1) + 1)
      sdstep[probs<lower1] <- sdstep[probs<lower1] * 1 / ( 2 - (probs[probs<lower1])/lower1 )
      # store sdstep
      sdstepseq[sdsteptel,1] <- sdstep
      sdsteptel <- sdsteptel + 1
    }

    store[i,] <- c(mu)

    #print(c(mu,tau2,sdstep))

  }

  return(list(postdraws = as.matrix(store[-(1:burnin),]), sdstepseq=sdstepseq))

}

# marginal likelihood computation for #1
log_marg_like_full <- function(yi, vi, tau2_min,
                               prior.mu,
                               prior.tau2,
                               start_mu, start_tau2, sdstep.mu, sdstep.tau2,
                               burnin, iters1, iters2){

  #estimate unconstrained marema model with proper prior for mu
  unc.draws.marglike <- gibbs_unc_prior(yi, vi, tau2_min,
                                        prior.mu,
                                        prior.tau2,
                                        start_mu, start_tau2, sdstep.mu, sdstep.tau2,
                                        burnin, iters1)

  #approximate posterior tau2 with a log normal for constructing proposal distribution
  para.tau2 <- fitdistr(unc.draws.marglike[[1]][,2]-tau2_min, "log-normal", lower = 0.0000001)$estimate
  para.tau2.IS <- c(para.tau2[1],para.tau2[2]*2)
  draws.IS.tau2 <- rlnorm(iters2, meanlog = para.tau2.IS[1], sdlog = para.tau2.IS[2]) + tau2_min

  #approximate posterior mu with a normal
  para.mu <- c(mean(unc.draws.marglike[[1]][,1]),sd(unc.draws.marglike[[1]][,1]))
  para.mu.IS <- para.mu
  para.mu.IS[2] <- para.mu.IS[2] * 2
  draws.IS.mu <- rnorm(iters2, mean = para.mu.IS[1], sd = para.mu.IS[2])

  #evaluate integrand of marginal likelihood for proposal draws
  logintegrands <- unlist(lapply(1:iters2,function(s){
    sum(dnorm(yi,mean=draws.IS.mu[s],sd=sqrt(vi+draws.IS.tau2[s]),log=TRUE)) +
      prior.mu(draws.IS.mu[s],draws.IS.tau2[s],log=TRUE) + prior.tau2(draws.IS.tau2[s],log=TRUE) -
      dnorm(draws.IS.mu[s], mean = para.mu.IS[1], sd = para.mu.IS[2],log=TRUE) -
      dlnorm(draws.IS.tau2[s]-tau2_min, meanlog = para.tau2.IS[1], sdlog = para.tau2.IS[2], log=TRUE)
  }))
  logintegrands <- logintegrands[!is.na(logintegrands)]

  logmuu <- log(mean(exp(logintegrands-max(logintegrands)))) + max(logintegrands)

  return(list(log.marg.like = logmuu, post.draws = unc.draws.marglike[[1]], sdstepseq = unc.draws.marglike$sdstepseq))

}

# marginal likelihood computation for #2
log_marg_like_cond.mu <- function(yi, vi, tau2_min,
                                  muIN, prior.tau2,
                                  start_tau2, sdstep.tau2, burnin, iters1, iters2){

  draws.marglike.cond.mu <- gibbs_cond.mu_prior(yi, vi, tau2_min,
                                                muIN, prior.tau2=prior.tau2,
                                                start_tau2, sdstep.tau2, burnin, iters1)

  #approximate posterior tau2 with a shifted log-normal as proposal distribution
  para.tau2 <- fitdistr(draws.marglike.cond.mu[[1]][,1]-tau2_min, "log-normal", lower = 0.000001)$estimate
  para.tau2.IS <- c(para.tau2[1],para.tau2[2]*2)
  draws.IS.tau2 <- rlnorm(iters2, meanlog = para.tau2.IS[1], sdlog = para.tau2.IS[2]) + tau2_min

  #evaluate integrand of marginal likelihood for proposal draws
  logintegrands <- unlist(lapply(1:iters2,function(s){
    sum(dnorm(yi,mean=muIN,sd=sqrt(vi+draws.IS.tau2[s]),log=TRUE)) +
      prior.tau2(draws.IS.tau2[s],log=TRUE) -
      dlnorm(draws.IS.tau2[s]-tau2_min, meanlog = para.tau2.IS[1], sdlog = para.tau2.IS[2], log=TRUE)
  }))
  logintegrands <- logintegrands[!is.na(logintegrands)]
  logmuu_cond.mu <- log(mean(exp(logintegrands-max(logintegrands)))) + max(logintegrands)

  return(list(log.marg.like = logmuu_cond.mu, post.draws = draws.marglike.cond.mu[[1]], sdstepseq = draws.marglike.cond.mu[[2]]))

}

# marginal likelihood computation for #3
log_marg_like_cond.tau2 <- function(yi, vi, tau2IN,
                                    prior.mu,
                                    start_mu, sdstep.mu, burnin, iters1, iters2){

  draws.marglike.cond.tau2 <- gibbs_cond.tau_prior(yi, vi, tau2IN,
                                                   prior.mu, start_mu, sdstep.mu,
                                                   burnin, iters1)

  #approximate posterior mu with a normal
  para.mu <- c(mean(draws.marglike.cond.tau2[[1]][,1]),sd(draws.marglike.cond.tau2[[1]][,1]))
  para.mu.IS <- para.mu
  para.mu.IS[2] <- para.mu.IS[2] * 3
  draws.IS.mu <- rnorm(iters2, mean = para.mu.IS[1], sd = para.mu.IS[2])

  #evaluate integrand of marginal likelihood for proposal draws
  logintegrands <- unlist(lapply(1:iters2,function(s){
    sum(dnorm(yi,mean=draws.IS.mu[s],sd=sqrt(vi+tau2IN),log=TRUE)) +
      prior.mu(draws.IS.mu[s],tau2IN,log=TRUE) -
      dnorm(draws.IS.mu[s], mean = para.mu.IS[1], sd = para.mu.IS[2],log=TRUE)
  }))
  logintegrands <- logintegrands[!is.na(logintegrands)]

  logmu0 <- log(mean(exp(logintegrands-max(logintegrands)))) + max(logintegrands)

  return(list(log.marg.like = logmu0, post.draws = draws.marglike.cond.tau2[[1]], sdstepseq = draws.marglike.cond.tau2[[2]]))

}


# test different priors for tau
BF_rma.uni <- function(x,
                       hypothesis = NULL,
                       prior.hyp.explo = NULL,
                       prior.hyp.conf = NULL,
                       prior.hyp = NULL,
                       complement = TRUE,
                       log = FALSE,
                       BF.type,
                       prior.tau2 = 1,
                       iter = 2e4,
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

  if(!exists("BF.type")){
    stop("The argument 'BF.type' is missing. See documentation. See ?BF")
  }

  logIN <- log

  # check proper usage of argument 'prior.hyp.conf' and 'prior.hyp.explo'
  if(!is.null(prior.hyp.conf)){
    prior.hyp <- prior.hyp.conf
  }
  prior.hyp.explo <- process.prior.hyp.explo(prior_hyp_explo = prior.hyp.explo, model = x)

  ### Extract effect sizes and sampling variances from the metafor object
  yi <- c(x$yi)
  vi <- x$vi
  ni <- x$ni #only used in case of unit-information prior: BF.type = "unit.info"
  N <- sum(ni)
  K <- length(vi)

  wi <- 1/vi # Weights equal-effects model
  typ_vi <- sum(wi*(length(wi)-1))/(sum(wi)^2 - sum(wi^2)) # Typical within-study sampling variance
  est.mu <- sum(yi*wi) / sum(wi)
  se.mu <- x$se

  est.tau2 <- x$tau2
  se.tau2 <- x$se.tau2 # check if we can use it as step sd for posterior sampling of tau
  if(is.na(se.tau2)){se.tau2 <- .5}
  tau2.min <- -min(vi)

  # prior tau2
  if(prior.tau2 == 1){
    prior.tau2.Jeffreys1 <- function(tau2.arg,log=FALSE){
      ptau2 <- -1/K*sum(log(vi+tau2.arg))
      #ptau2 <- .5*log(sum(1/(vi+tau2.arg)^2)) # Tibshirani (1989)
      #ptau2 <- 0 #uniform prior
      ifelse(log,ptau2,exp(ptau2))
    }
  }else if(prior.tau2 == 2){
    prior.tau2.Jeffreys1 <- function(tau2.arg,log=FALSE){
      #ptau2 <- -1/K*sum(log(vi+tau2.arg))
      ptau2 <- .5*log(sum(1/(vi+tau2.arg)^2)) # Tibshirani (1989)
      #ptau2 <- 0 #uniform prior
      ifelse(log,ptau2,exp(ptau2))
    }
  }else if(prior.tau2 == 3){
    prior.tau2.Jeffreys1 <- function(tau2.arg,log=FALSE){
      #ptau2 <- -1/K*sum(log(vi+tau2.arg))
      #ptau2 <- .5*log(sum(1/(vi+tau2.arg)^2)) # Tibshirani (1989)
      ptau2 <- 0 #uniform prior
      ifelse(log,ptau2,exp(ptau2))
    }
  }else{
    prior.tau2.Jeffreys1 <- function(tau2.arg,log=FALSE){
      #ptau2 <- -1/K*sum(log(vi+tau2.arg))
      #ptau2 <- .5*log(sum(1/(vi+tau2.arg)^2)) # Tibshirani (1989)
      #ptau2 <- 0 #uniform prior
      #ptau2 <- dnorm(tau2.arg,mean=tau2.min,sd=1e3,log=log) * 2
      ptau2 <- tau2.arg - tau2.min
      ifelse(log,ptau2,exp(ptau2))
    }
  }

  if(is(BF.type,"character")){
    if(BF.type == "stand.effect"){
      # conjugate normal prior assuming average effects of about .5
      prior.mu <- function(x,tau2,log=FALSE){dnorm(x, mean = 0, sd = .5, log = log)}
      prior.muGR0 <- .5
      bayesfactor.name <- "Bayes factor based on a normal prior (mu ~ norm(mean=0, sd=0.5))"
      parameter.name <- "standardized effect"
    }else if(BF.type == "log.odds"){
      # Student t prior which approximates the implied distribution of the log odds ratio based on uniform success probabilities
      prior.mu <- function(x,tau2,log=FALSE){
        dt1 <- dt(x/2.36,df=13.1,log=log)-log(2.36)
        ifelse(log,dt1,exp(dt1))
      }
      prior.muGR0 <- .5
      bayesfactor.name <- "Bayes factor based on uniform priors for success probabilities (log.odds ~ t(0,2.36,13.1))."
      parameter.name <- "log odds"
    }else if(BF.type == "correlation"){
      # logistic prior for the Fisher transformed correlation corresponding to a uniform prior for the correlation in (-1,1)
      prior.mu <- function(x,tau2,log=FALSE){dlogis(x, scale = .5, log = log)}
      prior.muGR0 <- .5
      bayesfactor.name <- "Bayes factor based on a uniform prior for the correlation in (-1,1) (Fisher(cor)~logis(0.5))."
      parameter.name <- "correlation"
    }else if(BF.type == "unit.info"){
      if(N==0){
        stop("In order to use a unit-information prior, the sample sizes need to be part of the element 'ni'
              of the 'rma.uni' object (e.g., using the 'ni' argument of the 'rma' function).")
      }
      prior.mu <- function(x,tau2,log=FALSE){dnorm(x, mean = 0, sd = sqrt(N/sum(1/(vi+tau2))), log = log)}
      prior.muGR0 <- .5
      bayesfactor.name <- "Bayes factor based on unit information prior."
      parameter.name <- "general parameter"
    }else{
      stop("The argument 'BF.type' is not correctly specified for an object of type 'rma.uni'. See documentation. ?BF")
    }
  }else if(is(BF.type,"prior")){
    prior.mu <- function(x,tau2,log=FALSE){BF.type(x,log)}
    seq1 <- seq(0,1e3,length=1e5)
    seq1 <- seq1+(seq1[2]-seq1[1])/2
    prior.muGR0 <- sum(prior.mu(seq1,0))*(seq1[2]-seq1[1]) #riemann estimate of prior prob that mu > 0.
    bayesfactor.name <- paste0("Bayes factor based on manually chosen prior")
    parameter.name <- "general parameter"
    message("Be sure that the specified prior is not truncated in a specific interval")
  }else{
    stop("The argument 'BF.type' is not correctly specified for an object of type 'rma.uni'. See documentation. ?BF")
  }

  #exploratory testing
  if (x$method != "EE" & x$method != "FE") ### random effects / marema model
  {

    #### MAREMA MODEL

    # initial values for MCMC run
    start_mu <- mean(yi)
    sdstep.mu <- x$se
    start_tau2 <- ifelse(x$tau2>0, x$tau2,.05)
    sdstep.tau2 <- x$se.tau2

    # marginal likelihood full unconstrained model. marema.unc.post.draws
    marg.like.marema.full <- log_marg_like_full(yi, vi, tau2_min=tau2.min,
                                                prior.mu=prior.mu,
                                                prior.tau2=prior.tau2.Jeffreys1,
                                                start_mu=start_mu, start_tau2=start_tau2,
                                                sdstep.mu=sdstep.mu, sdstep.tau2=sdstep.tau2,
                                                burnin = round(iter/2), iters1=iter, iters2=iter)

    # marginal likelihood under mu = 0
    marg.like.marema.muEQ0 <- log_marg_like_cond.mu(yi, vi, tau2_min=tau2.min,
                                                    muIN=0, prior.tau2=prior.tau2.Jeffreys1,
                                                    start_tau2=start_tau2+start_mu^2, sdstep.tau2=sdstep.tau2,
                                                    burnin=iter, iters1=iter, iters2=iter)

    # marginal likelihoods when constraining mu
    # mu = 0
    marg.like.marema.muEQ0 <- marg.like.marema.muEQ0[[1]]
    # mu > 0
    marg.like.marema.muGR0 <- marg.like.marema.full[[1]] + log(mean(marg.like.marema.full$post.draws[,1]>0)/prior.muGR0)
    # mu < 0
    marg.like.marema.muSM0 <- marg.like.marema.full[[1]] + log(mean(marg.like.marema.full$post.draws[,1]<0)/(1-prior.muGR0))

    #### random effects MODEL

    # marginal likelihood full unconstrained model
    marg.like.ranef.full <- log_marg_like_full(yi, vi, tau2_min=0,
                                               prior.mu = prior.mu,
                                               prior.tau2 = prior.tau2.Jeffreys1,
                                               start_mu = start_mu, start_tau2 = ifelse(start_tau2<0,sdstep.tau2,start_tau2),
                                               sdstep.mu = sdstep.mu, sdstep.tau2 = sdstep.tau2,
                                               burnin=round(iter/2), iters1=iter, iters2=iter)

    # marginal likelihood under mu = 0
    marg.like.ranef.muEQ0 <- log_marg_like_cond.mu(yi, vi, tau2_min=0,
                                                   muIN=0, prior.tau2=prior.tau2.Jeffreys1,
                                                   start_tau2=ifelse(start_tau2<0,sdstep.tau2,start_tau2)+start_mu^2, sdstep.tau2=sdstep.tau2,
                                                   burnin = round(iter/2), iters1=iter, iters2=iter)

    marg.like.ranef.muEQ0 <- marg.like.ranef.muEQ0[[1]]
    # mu > 0
    marg.like.ranef.muGR0 <- marg.like.ranef.full[[1]] + log(mean(marg.like.ranef.full$post.draws[,1]>0)/prior.muGR0)
    # mu < 0
    marg.like.ranef.muSM0 <- marg.like.ranef.full[[1]] + log(mean(marg.like.ranef.full$post.draws[,1]<0)/(1-prior.muGR0))

    ### Compute Bayes factors model vs. unconstrained mu and I^2 (marema)
    BF0u_ranef <- marg.like.ranef.muEQ0 - marg.like.ranef.full[[1]]
    BF1u_ranef <- marg.like.ranef.muSM0 - marg.like.ranef.full[[1]]
    BF2u_ranef <- marg.like.ranef.muGR0 - marg.like.ranef.full[[1]]
    BF0u_marema <- marg.like.marema.muEQ0 - marg.like.marema.full[[1]]
    BF1u_marema <- marg.like.marema.muSM0 - marg.like.marema.full[[1]]
    BF2u_marema <- marg.like.marema.muGR0 - marg.like.marema.full[[1]]

    BFtu_exploratory_marema <- matrix(c(BF0u_marema,BF1u_marema,BF2u_marema),nrow=1)
    rowmax <- apply(BFtu_exploratory_marema,1,max)
    norm_BF_explo <- exp(BFtu_exploratory_marema - rowmax %*% t(rep(1,ncol(BFtu_exploratory_marema)))) *
      (rep(1,nrow(BFtu_exploratory_marema)) %*% t(prior.hyp.explo[[1]]))
    PHP_exploratory_marema <- norm_BF_explo / apply(norm_BF_explo,1,sum)

    BFtu_exploratory_ranef <- matrix(c(BF0u_ranef,BF1u_ranef,BF2u_ranef),nrow=1)
    rowmax <- apply(BFtu_exploratory_ranef,1,max)
    norm_BF_explo <- exp(BFtu_exploratory_ranef - rowmax %*% t(rep(1,ncol(BFtu_exploratory_ranef)))) *
      (rep(1,nrow(BFtu_exploratory_ranef)) %*% t(prior.hyp.explo[[1]]))
    PHP_exploratory_ranef <- norm_BF_explo / apply(norm_BF_explo,1,sum)

    BFtu_exploratory <- rbind(BFtu_exploratory_marema,BFtu_exploratory_ranef)
    PHP_exploratory <- rbind(PHP_exploratory_marema,PHP_exploratory_ranef)

    colnames(BFtu_exploratory) <- c("=0","<0",">0")
    colnames(PHP_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
    row.names(BFtu_exploratory) <- row.names(PHP_exploratory) <- c("mu (marema)", "mu (ranef)")

    if(logIN == FALSE){
      BFtu_exploratory <- exp(BFtu_exploratory)
    }

    #store descriptives of estimates
    estimates_marema <- cbind(apply(marg.like.marema.full$post.draws,2,mean),
                              apply(marg.like.marema.full$post.draws,2,median),
                              apply(marg.like.marema.full$post.draws,2,quantile,.025),
                              apply(marg.like.marema.full$post.draws,2,quantile,.975),
                              apply(marg.like.marema.full$post.draws,2,function(x){mean(x>0)}))
    row.names(estimates_marema) <- c("mu      (marema)", "tau2    (marema)")
    estimates_ranef <- cbind(apply(marg.like.ranef.full$post.draws,2,mean),
                             apply(marg.like.ranef.full$post.draws,2,median),
                             apply(marg.like.ranef.full$post.draws,2,quantile,.025),
                             apply(marg.like.ranef.full$post.draws,2,quantile,.975),
                             apply(marg.like.ranef.full$post.draws,2,function(x){mean(x>0)}))
    row.names(estimates_ranef) <- c("mu      (ranef)", "tau2    (ranef)")
    colnames(estimates_marema) <- colnames(estimates_ranef) <- c("mean","median","2.5%","97.5%","Pr(>0)")
    uncestimates <- estimates_marema
    # estimates of separate group means
    #under marema model
    marema_post.draws_mu <- marg.like.marema.full$post.draws[,1]
    marema_post.draws_tau2 <- marg.like.marema.full$post.draws[,2]
    tau2draws_trunc <- marg.like.marema.full$post.draws[marema_post.draws_tau2>0,2]
    mudraws_trunc <- marg.like.marema.full$post.draws[marema_post.draws_tau2>0,1]
    numzerotau2 <- sum(marema_post.draws_tau2 <= 0)
    mudraws_studies <- do.call(cbind,lapply(1:length(vi),function(s){
      var_s <- 1/(1/vi[s] + 1/tau2draws_trunc)
      mean_s <- (yi[s]/vi[s] + mudraws_trunc/tau2draws_trunc) * var_s
      draws_s <- c(rnorm(length(mudraws_trunc),mean=mean_s,sd=sqrt(var_s)),
                   marema_post.draws_mu[marema_post.draws_tau2<0])
      c(draws_s,mean(draws_s),median(draws_s),quantile(draws_s,.025),quantile(draws_s,.975),mean(draws_s>0))
    }))
    uncestimates <- rbind(uncestimates,t(mudraws_studies[length(marema_post.draws_mu)+1:5,]))
    row.names(uncestimates)[-(1:2)] <- paste0("theta_",1:ncol(mudraws_studies)," (marema)")
    #under random effects model
    uncestimates <- rbind(uncestimates,estimates_ranef)
    ranef_post.draws_mu <- marg.like.ranef.full$post.draws[,1]
    ranef_post.draws_tau2 <- marg.like.ranef.full$post.draws[,2]
    mudraws_studies <- do.call(cbind,lapply(1:length(vi),function(s){
      var_s <- 1/(1/vi[s] + 1/ranef_post.draws_tau2)
      mean_s <- (yi[s]/vi[s] + ranef_post.draws_mu/ranef_post.draws_tau2) * var_s
      draws_s <- rnorm(length(ranef_post.draws_mu),mean=mean_s,sd=sqrt(var_s))
      c(draws_s,mean(draws_s),median(draws_s),quantile(draws_s,.025),quantile(draws_s,.975),mean(draws_s>0))
    }))
    uncestimates <- rbind(uncestimates,t(mudraws_studies[length(ranef_post.draws_mu)+1:5,]))
    row.names(uncestimates)[-(1:(4+length(vi)))] <- paste0("theta_",1:ncol(mudraws_studies)," (ranef)")

    #store posterior draws based on noninformative improper prior
    post.draws <- list()
    post.draws$marema <- marg.like.marema.full$post.draws
    post.draws$ranef <- marg.like.ranef.full$post.draws

  } else if (x$method == "EE" | x$method == "FE") { ### common effect model

    # testing
    marg.like.comef.muEQ0 <- log_marg_like_cond.tau2(yi, vi, tau2IN=0,
                                                     prior.mu,
                                                     start_mu=est.mu, sdstep.mu=se.mu,
                                                     burnin = round(iter/2), iters1 = iter, iters2 = iter)

    ### log marginal likelihoods
    # mu = 0
    marg.like.muEQ0 <- sum(dnorm(yi,mean=0,sd=sqrt(vi+0),log=TRUE))
    # mu > 0
    marg.like.muGR0 <- marg.like.comef.muEQ0[[1]] + log(mean(marg.like.comef.muEQ0[[2]]>0)/prior.muGR0)
    # mu < 0
    marg.like.muSM0 <- marg.like.comef.muEQ0[[1]] + log(mean(marg.like.comef.muEQ0[[2]]<0)/(1-prior.muGR0))

    ### Compute Bayes factors model vs. unconstrained model
    BF0u_comef <- marg.like.muEQ0 - marg.like.comef.muEQ0[[1]]
    BF1u_comef <- marg.like.muSM0 - marg.like.comef.muEQ0[[1]]
    BF2u_comef <- marg.like.muGR0 - marg.like.comef.muEQ0[[1]]

    BFtu_exploratory <- matrix(c(BF0u_comef,BF1u_comef,BF2u_comef),nrow=1)
    rowmax <- apply(BFtu_exploratory,1,max)
    norm_BF_explo <- exp(BFtu_exploratory - rowmax %*% t(rep(1,ncol(BFtu_exploratory)))) *
      (rep(1,nrow(BFtu_exploratory)) %*% t(prior.hyp.explo[[1]]))
    PHP_exploratory <- norm_BF_explo / apply(norm_BF_explo,1,sum)
    colnames(BFtu_exploratory) <- c("=0","<0",">0")
    colnames(PHP_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
    row.names(BFtu_exploratory) <- row.names(PHP_exploratory) <- c("mu (comef)")

    #descriptives of estimates
    uncestimates <- cbind(mean(marg.like.comef.muEQ0$post.draws), median(marg.like.comef.muEQ0$post.draws),
                          quantile(marg.like.comef.muEQ0$post.draws,.025),quantile(marg.like.comef.muEQ0$post.draws,.975),
                          mean(marg.like.comef.muEQ0$post.draws>0))
    colnames(uncestimates) <- c("mean","median","2.5%","97.5%","Pr(>0)")
    row.names(uncestimates) <- "mu"

    #store posterior draws based on noninformative improper prior
    post.draws <- marg.like.comef.muEQ0$post.draws
    colnames(post.draws) <- row.names(uncestimates)

    if(logIN == FALSE){
      BFtu_exploratory <- exp(BFtu_exploratory)
    }

  }

  if(!is.null(hypothesis)){
    message("The 'hypothesis' argument for this model class is not yet supported. Please inform us about the types of hypotheses you wish to test :)")
  }

  BFtu_confirmatory <- PHP_confirmatory <- BFmatrix_confirmatory <- BFtable <-
    priorprobs <- hypotheses <- NULL

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
    bayesfactor=bayesfactor.name,
    parameter=parameter.name,
    log=logIN,
    call=match.call(),
    post.draws=post.draws
  )

  class(BF_out) <- "BF"

  BF_out

}





