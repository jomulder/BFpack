library(BFpack)
# example analysis
mtcars0 <- mtcars
mtcars0$vs[1:6] <- 2
mtcars0$vs <- as.factor(mtcars0$vs)
mtcars0$am <- as.factor(mtcars0$am)
# standardize nonfactor variables
mtcars0[,(names(mtcars0)!="vs")*(names(mtcars0)!="am")==1] <-
  scale(mtcars0[,(names(mtcars0)!="vs")*(names(mtcars0)!="am")==1])

# univariate regression
lm1 <- lm(wt ~ -1 + disp + vs + hp + drat, mtcars)
constraints <- "disp > drat > hp ; hp > disp = 0"
BFreg1 <- BF(lm1)
BFreg1
BFreg1 <- BF(lm1, hypothesis = constraints)
BFreg1
# update for later
# BFreg2 <- BFupdate.lm(BFreg1,lm1)

# ANOVA
lm1 <- lm(wt ~ -1 + vs, mtcars0)
constraints <- "vs0 > vs1 > vs2 ; vs0 = vs1 = vs2"
BFreg1 <- BF(lm1,hypothesis=constraints,priorprob="default")
# BFreg2 <- BFupdate.lm(BFreg1,lm1)

# ANCOVA
lm1 <- lm(wt ~ -1 + disp + vs + hp + drat, mtcars0)
constraints <- "vs0 > vs1 > vs2 ; vs0 = vs1 = vs2"
BFreg1 <- BF(lm1,hypothesis=constraints,priorprob="default")
# BFreg2 <- BFupdate.lm(BFreg1,lm1)

# multivariate regression
lm1 <- lm(cbind(wt,mpg) ~ -1 + disp + vs + hp + drat, mtcars0)
constraints <- "disp.wt > disp.mpg ; hp.wt > disp.wt = 0"
BFreg1 <- BF(lm1)
BFreg1 <- BF(lm1,hypothesis=constraints)
BFreg1
# BFreg2 <- BFupdate.lm(BFreg1,lm1)

# Student t test
ttest <- t.test(-3:10)
BF(ttest)
BF(ttest,"exploratory")
BF(ttest,"mu_minus_mu0<0")


#mixed outcomes
mtcars0$vs <- as.factor(mtcars0$vs)
lm1 <- lm(cbind(qsec,vs) ~ -1 + disp + wt + hp + drat, mtcars0)
lm1 <- lm(cbind(qsec,gear) ~ -1 + disp + wt + hp + drat, mtcars0)
lm1$terms
lm1$coefficients
lm1$model
attr(terms(lm1),"factors")
attr(terms(lm1),"dataClasses")
attr(terms(lm1),"term.labels")
attr(terms(lm1),"specials")
attr(terms(lm1),"variables")[2]


attributes(terms(lm1))

terms(lm1)$factors
attributes(lm1,"dataClasses")

# test correlation analysis
mtcars$vss <- as.factor(mtcars$vs)
lm1 <- lm(cbind(mpg,cyl,disp) ~ -1 + wt + vss, mtcars)
constraints="mpg_with_cyl_in_vss0>mpg_with_disp_in_vss1>cyl_with_disp_in_vss0;mpg_with_cyl_in_vss0>mpg_with_disp_in_vss1=0"


summary(lm1)
model.matrix(lm1)

coefficients(lm1)
BFcorr(lm1,constraints="exploratory")

library(polycor)
my_data <- mtcars[, c(1,3,4,5,6,7)]
res <- cor(my_data)
res <- hetcor(my_data)

install.packages("Hmisc")
library("Hmisc")
res2 <- rcorr(as.matrix(my_data))
res2$P





# make a factor of variable site
sesamesim$site <- as.factor(sesamesim$site)
# execute an analysis of variance using lm() which, due to the -1, returns
# estimates of the means per group
anov <- lm(postnumb~site-1,sesamesim)
# take a look at the estimated means and their names
coef(anov)
# set a seed value
set.seed(100)
# use the names to formulate and test hypotheses with bain
results <- bain(anov, "site1=site2=site3=site4=site5; site2>site5>site1>
                site3>site4")




library(lmhyp)
fit <- lm(mpg ~ disp + hp + wt, data = mtcars)
Hyp <- "wt > disp > hp > 0; wt > hp > disp > 0"
result <- test_hyp(fit, Hyp, mcrep = 1000000)
result
