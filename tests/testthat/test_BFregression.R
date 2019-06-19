# # test BF on lm object
#
# # library(BFpack)
# # # example analysis
# # mtcars$vsf <- as.factor(mtcars$vs)
# # mtcars$amf <- as.factor(mtcars$am)
# # mtcars$gearf <- as.factor(mtcars$gear)
# # # standardize nonfactor variables
# # # mtcars0[,(names(mtcars0)!="vs")*(names(mtcars0)!="am")==1] <-
# # #   scale(mtcars0[,(names(mtcars0)!="vs")*(names(mtcars0)!="am")==1])
# #
# # univariate regression
# lm1 <-  lm(wt ~ disp + drat + hp, mtcars)
# BF(lm1)
# BF(lm1,hypothesis="disp=drat=0;disp>drat>0")
#
# #multivariate regression
# lm1 <-  lm(cbind(wt,disp) ~ drat + hp, mtcars)
# BF(lm1)
# BF(lm1,hypothesis="drat_on_wt=drat_on_disp=0;drat_on_wt=drat_on_disp>0")
#
# # lm1 <- aov(wt ~ -1 + gearf + drat + disp + cyl+ vsf + amf:vsf:drat:vsf + carb:gearf + qsec + hp, mtcars)
# # BF(lm1)
# # lm1 <- aov(wt  ~ -1 + gearf + amf:vsf:disp, mtcars)
# # BF(lm1)
# #
# # constraints <- "(disp , drat) > (hp,vs)"
# # constraints <- "(disp , drat) > hp; disp = drat = hp; drat > hp"
# # constraints <- "(disp, drat) > hp; disp = drat = hp; drat > hp"
# # BF1 <- BF(x=lm1, hypothesis = constraints)
# # BF1$relative_complexity
# #
# # # update for later
# # # BFreg2 <- BFupdate.lm(BFreg1,lm1)
# #
# # # ANOVA
# # lm1 <- lm(wt ~ -1 + vs, mtcars0)
# # constraints <- "vs0 > vs1 > vs2 ; vs0 = vs1 = vs2"
# # BFreg1 <- BF(lm1,hypothesis=constraints,priorprob="default")
# # # BFreg2 <- BFupdate.lm(BFreg1,lm1)
# #
# # # ANCOVA
# # lm1 <- lm(wt ~ -1 + disp + vs + hp + drat, mtcars0)
# # constraints <- "vs0 > vs1 > vs2 ; vs0 = vs1 = vs2"
# # BFreg1 <- BF(lm1,hypothesis=constraints,priorprob="default")
# # # BFreg2 <- BFupdate.lm(BFreg1,lm1)
# #
# # # multivariate regression
# # lm1 <- lm(cbind(wt,mpg) ~ 1 + disp + vs + hp + drat, mtcars)
# # BF(lm1)
# # constraints <- "(hp_on_wt,disp_on_wt) > (drat_on_mpg,disp_on_mpg) ; hp_on_wt > disp_on_wt = 0"
# # BF(lm1,hypothesis=constraints)
# # BFreg1
# # # BFreg2 <- BFupdate.lm(BFreg1,lm1)
# #
# # # Student t test
# # ttest <- t.test(-3:10)
# # BF(ttest)
# # BF(ttest,"exploratory")
# # BF(ttest,"mu_minus_mu0<0")
# #
# # #mixed outcomes
# # mtcars0$vs <- as.factor(mtcars0$vs)
# # lm1 <- lm(cbind(qsec,vs) ~ -1 + disp + wt + hp + drat, mtcars0)
# # lm1 <- lm(cbind(qsec,gear) ~ -1 + disp + wt + hp + drat, mtcars0)
# # lm1$terms
# # lm1$coefficients
# # lm1$model
# # attr(terms(lm1),"factors")
# # attr(terms(lm1),"dataClasses")
# # attr(terms(lm1),"term.labels")
# # attr(terms(lm1),"specials")
# # attr(terms(lm1),"variables")[2]
# #
# #
# # attributes(terms(lm1))
# #
# # terms(lm1)$factors
# # attributes(lm1,"dataClasses")
# #
# # # test correlation analysis
# # mtcars$vss <- as.factor(mtcars$vs)
# # lm1 <- lm(cbind(mpg,cyl,disp) ~ 1 + wt, mtcars)
# # constraints="mpg_with_cyl>disp_with_mpg>cyl_with_disp;mpg_with_cyl>disp_with_mpg=0"
# # BF(x=lm1,parameter="correlation")
# #
# # constraints="cyl_with_mpg > disp_with_mpg > 0"
# #
# # constraints="mpg_with_cyl_in_vss0>mpg_with_disp_in_vss1>cyl_with_disp_in_vss0;mpg_with_cyl_in_vss0>mpg_with_disp_in_vss1=0"
# #
# # constraints="(wt_on_mpg, wt_on_cyl) > (vss1_on_mpg,vss1_on_cyl) & (Intercept)_on_mpg = wt_on_mpg"
# #
# # constraints="(wt_on_mpg, wt_on_cyl) > (vss1_on_mpg,vss1_on_cyl) & vss1_on_mpg = wt_on_mpg"
# #
# # BF(lm1)
# # BF(lm1,hypothesis=constraints)
# # BF(x=lm1,parameter="correlation")
# #
# # BFcorr(lm1,constraints)
# #
# #
# # summary(lm1)
# # model.matrix(lm1)
# #
# # coefficients(lm1)
# # BFcorr(lm1,constraints="exploratory")
# #
# #
# #
# #
# #
# #
# #
# #
# # # make a factor of variable site
# # sesamesim$site <- as.factor(sesamesim$site)
# # # execute an analysis of variance using lm() which, due to the -1, returns
# # # estimates of the means per group
# # anov <- lm(postnumb~site-1,sesamesim)
# # # take a look at the estimated means and their names
# # coef(anov)
# # # set a seed value
# # set.seed(100)
# # # use the names to formulate and test hypotheses with bain
# # results <- bain(anov, "site1=site2=site3=site4=site5; site2>site5>site1>
# #                 site3>site4")
# #
# #
# # library(lmhyp)
# # fit <- lm(mpg ~ disp + hp + wt, data = mtcars)
# # Hyp <- "wt > disp > hp > 0; wt > hp > disp > 0"
# # result <- test_hyp(fit, Hyp, mcrep = 1000000)
# # result
# #
# #
# # npk.aov <- aov(yield ~ -1 + block + N*P*K, npk)
# # BF(npk.aov)
# #
# # Xmat <- model.matrix(npk.aov)
# # solve(t(Xmat)%*%Xmat)
# #
# # op <- options(contrasts = c("contr.helmert", "contr.poly"))
# # npk.aov <- aov(yield ~ block + N*P*K, npk)
# # summary(npk.aov)
# # coefficients(npk.aov)
# #
# #
# #
# #
# # cor1 <- polycor::hetcor(mvtnorm::rmvnorm(30,mean=rep(0,4),sigma=diag(4)+1))
# # constraints <- c("X2_with_X1>X1_with_X3>X1_with_X4;X2_with_X1=X3_with_X1=X4_with_X1>0;
# #   X2_with_X1=X3_with_X1=X4_with_X1=0")
# # x <- cor1
# # BF(cor1,constraints)
# #
# #
# # # ?manova test
# # npk2 <- within(npk, foo <- rnorm(24))
# # npk2.aov <- manova(cbind(yield, foo) ~ N*K, npk2)
# # summary(npk2.aov)
# # BF1 <- BF(npk2.aov)
# # BF1 <- BF(npk2.aov,"N1_on_yield>K1_on_yield>0")
# # class(npk2.aov)
# #
# # npk2.aov$coefficients
# #
# # npk2.aov$coefficients
# # model.matrix(npk2.aov)
# #
# # x=npk2.aov
# #
# #
# #
# #
# # devtools::install_github("jomulder/BFpack")
# # library(BFpack)
# # fit <- lm(cbind(mpg,disp) ~ 1 + hp + wt, data = mtcars)
# # BF(fit)
# # BF(fit,"hp_on_mpg>hp_on_disp;hp_on_mpg=hp_on_disp")
#
#
#
#
