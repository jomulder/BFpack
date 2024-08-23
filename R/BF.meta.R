
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
  est.mu <- sum(yi*wi) / sum(wi)
  se.mu <- x$se

  #typ_vi <- sum(wi*(length(wi)-1))/(sum(wi)^2 - sum(wi^2)) # Typical within-study sampling variance
  est.tau2 <- x$tau2
  se.tau2 <- x$se.tau2 # check if we can use it as step sd for posterior sampling of tau
  if(is.na(se.tau2)){se.tau2 <- .5}
  tau2.min <- -min(vi)

  # prior tau2
  prior.tau2.Jeffreys1 <- function(tau2.arg){prod(1/(vi+tau2.arg)^(1/K))} # Berger & Deely (1988)
  #prior.tau2.Jeffreys2 <- function(tau2.arg){sqrt(sum(1/(vi+tau2.arg)^2))} # Tibshirani (1989)
  #prior.tau2.uniform <- function(tau2.arg){1} # Tibshirani (1989)

  if(is(BF.type,"character")){
    if(BF.type == "stand.effect"){
      # conjugate normal prior assuming average effects of about .5
      prior.mu <- function(x,tau2){dnorm(x, mean = 0, sd = .5)}
      prior.muGR0 <- .5
      bayesfactor.name <- "Bayes factor based on a normal prior (mu ~ norm(mean=0, sd=0.5))"
      parameter.name <- "standardized effect"
    }else if(BF.type == "log.odds"){
      # Student t prior which approximates the implied distribution of the log odds ratio based on uniform success probabilities
      prior.mu <- function(x,tau2){dt(x/2.36,df=13.1)/2.36}
      prior.muGR0 <- .5
      bayesfactor.name <- "Bayes factor based on uniform priors for success probabilities (log.odds ~ t(0,2.36,13.1))."
      parameter.name <- "log odds"
    }else if(BF.type == "correlation"){
      # logistic prior for the Fisher transformed correlation corresponding to a uniform prior for the correlation in (-1,1)
      prior.mu <- function(x,tau2){dlogis(x, scale = .5)}
      prior.muGR0 <- .5
      bayesfactor.name <- "Bayes factor based on a uniform prior for the correlation in (-1,1) (Fisher(cor)~logis(0.5))."
      parameter.name <- "correlation"
    }else if(BF.type == "unit.info"){
      if(N==0){
        stop("In order to use a unit-information prior, the sample sizes need to be part of the element 'ni'
              of the 'rma.uni' object (e.g., using the 'ni' argument of the 'rma' function).")
      }
      prior.mu <- function(x,tau2){dnorm(x, mean = 0, sd = sqrt(N/sum(1/(vi+tau2))))}
      prior.muGR0 <- .5
      bayesfactor.name <- "Bayes factor based on unit information prior."
      parameter.name <- "general parameter"
    }else{
      stop("The argument 'BF.type' is not correctly specified for an object of type 'rma.uni'. See documentation. ?BF")
    }
  }else if(is(BF.type,"prior")){
    prior.mu <- function(x,tau2){BF.type(x)}
    seq1 <- seq(0,1e3,length=1e5)
    seq1 <- seq1+(seq1[2]-seq1[1])/2
    prior.muGR0 <- sum(prior.mu(seq1,0))*(seq1[2]-seq1[1]) #riemann estimate of prior prob that mu > 0.
    bayesfactor.name <- paste0("Bayes factor based on manually chosen prior")
    parameter.name <- "general parameter"
    message("Currently the function only supports priors in [-Inf,Inf]. Truncated priors should not be used.")
  }else{
    stop("The argument 'BF.type' is not correctly specified for an object of type 'rma.uni'. See documentation. ?BF")
  }

  #exploratory testing
  if (x$method != "EE" & x$method != "FE") ### random effects / marema model
  {

    #### MAREMA MODEL

    # unconstrained draws for estimation based on flat prior for mu
    marema.unc.post.draws <- gibbs_unc_prior(yi=yi, vi=vi, tau2_min=tau2.min,
                                             prior.mu = function(mu.arg,tau2.arg){1},
                                             prior.tau2 = prior.tau2.Jeffreys1,
                                             start_mu=est.mu, start_tau2=est.tau2,
                                             sdstep.mu = se.mu, sdstep.tau2 = se.tau2,
                                             burnin = round(iter/2), iters = iter)
    #get posterior estimates for initial values of subsequent analysis
    start_mu <- mean(marema.unc.post.draws$postdraws[,1])
    sdstep.mu <- sd(marema.unc.post.draws$postdraws[,1])
    start_tau2 <- mean(marema.unc.post.draws$postdraws[,2])
    sdstep.tau2 <- sd(marema.unc.post.draws$postdraws[,2])

    # marginal likelihood full unconstrained model
    marg.like.marema.full <- log_marg_like_full(yi, vi, tau2_min=tau2.min,
                                                prior.mu=prior.mu,
                                                prior.tau2=prior.tau2.Jeffreys1,
                                                start_mu=start_mu, start_tau2=start_tau2,
                                                sdstep.mu=sdstep.mu, sdstep.tau2=sdstep.tau2,
                                                burnin = round(iter/2), iters1=iter, iters2=iter)

    # marginal likelihood under mu = 0
    marg.like.marema.muEQ0 <- log_marg_like_cond.mu(yi, vi, tau2_min=tau2.min,
                                                    muIN=0, prior.tau2=prior.tau2.Jeffreys1,
                                                    start_tau2=start_tau2, sdstep.tau2=sdstep.tau2,
                                                    burnin = round(iter/2), iters1=iter, iters2=iter)

    # marginal likelihoods when constraining mu
    # mu = 0
    marg.like.marema.muEQ0 <- marg.like.marema.muEQ0[[1]]
    # mu > 0
    marg.like.marema.muGR0 <- marg.like.marema.full[[1]] + log(mean(marg.like.marema.full$post.draws[,1]>0)/prior.muGR0)
    # mu < 0
    marg.like.marema.muSM0 <- marg.like.marema.full[[1]] + log(mean(marg.like.marema.full$post.draws[,1]<0)/(1-prior.muGR0))

    #### random effects MODEL

    # unconstrained draws for estimation based on flat prior for mu
    ranef.unc.post.draws <- gibbs_unc_prior(yi=yi, vi=vi, tau2_min=0,
                                            prior.mu = function(mu.arg,tau2.arg){1},
                                            prior.tau2 = prior.tau2.Jeffreys1,
                                            start_mu = start_mu, start_tau2 = ifelse(start_tau2<0,sdstep.tau2,start_tau2),
                                            sdstep.mu = sdstep.mu, sdstep.tau2 = sdstep.tau2,
                                            burnin = round(iter/2), iters = iter)

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
                                                   start_tau2=ifelse(start_tau2<0,sdstep.tau2,start_tau2), sdstep.tau2=sdstep.tau2,
                                                   burnin = round(iter/2), iters1=iter, iters2=iter)

    marg.like.ranef.muEQ0 <- marg.like.ranef.muEQ0[[1]]
    # mu > 0
    marg.like.ranef.muGR0 <- marg.like.ranef.full[[1]] + log(mean(marg.like.ranef.full$post.draws[,1]>0)/prior.muGR0)
    # mu < 0
    marg.like.ranef.muSM0 <- marg.like.ranef.full[[1]] + log(mean(marg.like.ranef.full$post.draws[,1]<0)/(1-prior.muGR0))

    ### Compute Bayes factors model vs. unconstrained mu and I^2 (marema)
    BF0u_ranef <- marg.like.ranef.muEQ0 - marg.like.ranef.full[[1]]
    BF1u_ranef <- marg.like.ranef.muGR0 - marg.like.ranef.full[[1]]
    BF2u_ranef <- marg.like.ranef.muSM0 - marg.like.ranef.full[[1]]
    BF0u_marema <- marg.like.marema.muEQ0 - marg.like.marema.full[[1]]
    BF1u_marema <- marg.like.marema.muGR0 - marg.like.marema.full[[1]]
    BF2u_marema <- marg.like.marema.muSM0 - marg.like.marema.full[[1]]

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
    estimates_marema <- cbind(apply(marema.unc.post.draws$postdraws,2,mean),
                              apply(marema.unc.post.draws$postdraws,2,median),
                              apply(marema.unc.post.draws$postdraws,2,quantile,.025),
                              apply(marema.unc.post.draws$postdraws,2,quantile,.975),
                              apply(marema.unc.post.draws$postdraws,2,function(x){mean(x>0)}))
    row.names(estimates_marema) <- c("mu      (marema)", "tau2    (marema)")
    estimates_ranef <- cbind(apply(ranef.unc.post.draws$postdraws,2,mean),
                             apply(ranef.unc.post.draws$postdraws,2,median),
                             apply(ranef.unc.post.draws$postdraws,2,quantile,.025),
                             apply(ranef.unc.post.draws$postdraws,2,quantile,.975),
                             apply(ranef.unc.post.draws$postdraws,2,function(x){mean(x>0)}))
    row.names(estimates_ranef) <- c("mu      (ranef)", "tau2    (ranef)")
    colnames(estimates_marema) <- colnames(estimates_ranef) <- c("mean","median","2.5%","97.5%","Pr(>0)")
    uncestimates <- estimates_marema
    # estimates of separate group means
    #under marema model
    marema_post.draws_mu <- marema.unc.post.draws$postdraws[,1]
    marema_post.draws_tau2 <- marema.unc.post.draws$postdraws[,2]
    tau2draws_trunc <- marema.unc.post.draws$postdraws[marema_post.draws_tau2>0,2]
    mudraws_trunc <- marema.unc.post.draws$postdraws[marema_post.draws_tau2>0,1]
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
    ranef_post.draws_mu <- ranef.unc.post.draws$postdraws[,1]
    ranef_post.draws_tau2 <- ranef.unc.post.draws$postdraws[,2]
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
    post.draws$marema <- marema.unc.post.draws$postdraws
    post.draws$ranef <- ranef.unc.post.draws$postdraws

  } else if (x$method == "EE" | x$method == "FE") { ### common effect model

    # estimation
    comef.unc.post.draws <- gibbs_cond.tau_prior(yi, vi, tau2IN = 0,
                                                 prior.mu,
                                                 start_mu=est.mu, sdstep.mu=se.mu,
                                                 burnin = round(iter/2), iters = iter)

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
    BF1u_comef <- marg.like.muGR0 - marg.like.comef.muEQ0[[1]]
    BF2u_comef <- marg.like.muSM0 - marg.like.comef.muEQ0[[1]]

    BFtu_exploratory <- matrix(c(BF0u_comef,BF1u_comef,BF2u_comef),nrow=1)
    rowmax <- apply(BFtu_exploratory,1,max)
    norm_BF_explo <- exp(BFtu_exploratory - rowmax %*% t(rep(1,ncol(BFtu_exploratory)))) *
      (rep(1,nrow(BFtu_exploratory)) %*% t(prior.hyp.explo[[1]]))
    PHP_exploratory <- norm_BF_explo / apply(norm_BF_explo,1,sum)
    colnames(BFtu_exploratory) <- c("=0","<0",">0")
    colnames(PHP_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
    row.names(BFtu_exploratory) <- row.names(PHP_exploratory) <- c("mu (comef)")

    #descriptives of estimates
    uncestimates <- cbind(mean(comef.unc.post.draws$postdraws), median(comef.unc.post.draws$postdraws),
                          quantile(comef.unc.post.draws$postdraws,.025),quantile(comef.unc.post.draws$postdraws,.975),
                          mean(comef.unc.post.draws$postdraws>0))
    colnames(uncestimates) <- c("mean","median","2.5%","97.5%","Pr(>0)")
    row.names(uncestimates) <- "mu"

    #store posterior draws based on noninformative improper prior
    post.draws <- comef.unc.post.draws$postdraws
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
    #BF_delta = BF_delta,
    #BF_rho = BF_rho,
    #PP_rho = PP_rho
  )

  class(BF_out) <- "BF"

  BF_out

}

#1. unconstrained sampling under marema using 'prior.mu' as prior for mu and a noninformative prior for tau in (prior.tau2,Inf)
gibbs_unc_prior <- function(yi, vi, tau2_min,
                            prior.mu, prior.tau2,
                            start_mu, start_tau2, sdstep.mu, sdstep.tau2, burnin, iters){

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
  sdstepseq <- matrix(0, nrow = (iters + burnin) / check1,2)
  sdstep <- c(sdstep.mu,sdstep.tau2)
  store <- matrix(NA, ncol=2, nrow = iters + burnin)
  colnames(store) <- colnames(acceptMat) <- colnames(sdstepseq) <- c("mu","tau2")

  for(i in 1:(burnin+iters)){
    #draw rho
    #draw candidate from truncated normal
    tau2_star <- rtnorm(1, mean = tau2, sd = sdstep[2], a = tau2_min)

    #evaluate Metropolis-Hastings acceptance probability
    can_vars <- vi+tau2_star
    # norm_cur <- 1 - pnorm(tau2_min,mean=tau2,sd=sdstep[2])
    # norm_can <- 1 - pnorm(tau2_min,mean=tau2_star,sd=sdstep[2])
    # norm_cur.div.norm_can <- ifelse(norm_cur == 0| norm_can == 0,1,norm_cur/norm_can)
    R_MH <- exp( sum(dnorm(yi,mean=mu,sd=sqrt(can_vars),log=TRUE)) -
                   sum(dnorm(yi,mean=mu,sd=sqrt(cur_vars),log=TRUE)) ) *
      (1 - pnorm(tau2_min,mean=tau2,sd=sdstep[2])) / (1 - pnorm(tau2_min,mean=tau2_star,sd=sdstep[2])) *
      prior.mu(mu,tau2_star) / prior.mu(mu,tau2) *
      prior.tau2(tau2_star) / prior.tau2(tau2)
    tau2 <- ifelse(runif(1) < R_MH, tau2_star, tau2)
    acceptMat[i,2] <- tau2_star == tau2
    cur_vars <- vi+tau2

    #draw mu
    mu_star <- rnorm(1, mean = mu, sd = sdstep[1])

    R_MH <- exp( sum(dnorm(yi,mean=mu_star,sd=sqrt(cur_vars),log=TRUE)) -
                   sum(dnorm(yi,mean=mu,sd=sqrt(cur_vars),log=TRUE)) ) *
      prior.mu(mu_star,tau2) / prior.mu(mu,tau2)
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

    #print(c(mu,tau2,sdstep))

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
    # norm_cur <- 1 - pnorm(tau2_min,mean=tau2,sd=sdstep[2])
    # norm_can <- 1 - pnorm(tau2_min,mean=tau2_star,sd=sdstep[2])
    # norm_cur.div.norm_can <- ifelse(norm_cur == 0| norm_can == 0,1,norm_cur/norm_can)
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

    #print(c(mu,tau2,sdstep))

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
      log(prior.mu(draws.IS.mu[s],draws.IS.tau2[s])) + log(prior.tau2(draws.IS.tau2[s])) -
      dnorm(draws.IS.mu[s], mean = para.mu.IS[1], sd = para.mu.IS[2],log=TRUE) -
      dlnorm(draws.IS.tau2[s]-tau2_min, meanlog = para.tau2.IS[1], sdlog = para.tau2.IS[2], log=TRUE)
  }))

  logmuu <- log(mean(exp(logintegrands-max(logintegrands)))) + max(logintegrands)
  logmuu

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
      log(prior.tau2(draws.IS.tau2[s])) -
      dlnorm(draws.IS.tau2[s]-tau2_min, meanlog = para.tau2.IS[1], sdlog = para.tau2.IS[2], log=TRUE)
  }))
  logmuu_cond.mu <- log(mean(exp(logintegrands-max(logintegrands)))) + max(logintegrands)
  logmuu_cond.mu

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
      log(prior.mu(draws.IS.mu[s],tau2IN)) -
      dnorm(draws.IS.mu[s], mean = para.mu.IS[1], sd = para.mu.IS[2],log=TRUE)
  }))

  logmu0 <- log(mean(exp(logintegrands-max(logintegrands)))) + max(logintegrands)

  return(list(log.marg.like = logmu0, post.draws = draws.marglike.cond.tau2[[1]], sdstepseq = draws.marglike.cond.tau2[[2]]))

}




