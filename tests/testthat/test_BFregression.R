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
lm1 <- lm(wt ~ -1 + disp + vs + hp + drat, mtcars0)
constraints <- "disp > drat > hp ; hp > disp = 0"
BFreg1 <- BF(lm1, hypothesis = constraints, prior="default")

BFreg2 <- BFregUpdate(BFreg1,lm1)

# ANOVA
lm1 <- lm(wt ~ -1 + vs, mtcars0)
constraints <- "vs0 > vs1 > vs2 ; vs0 = vs1 = vs2"
BFreg1 <- BFreg(lm1,constraints=constraints,priorprob="default")
BFreg2 <- BFregUpdate(BFreg1,lm1)

# ANCOVA
lm1 <- lm(wt ~ -1 + disp + vs + hp + drat, mtcars0)
constraints <- "vs0 > vs1 > vs2 ; vs0 = vs1 = vs2"
BFreg1 <- BFreg(lm1,constraints=constraints,priorprob="default")
BFreg2 <- BFregUpdate(BFreg1,lm1)

# multivariate regression
lm1 <- lm(cbind(wt,mpg) ~ -1 + disp + vs + hp + drat, mtcars0)
constraints <- "disp.wt > disp.mpg ; hp.wt > disp.wt = 0"
BFreg1 <- BFreg(lm1,constraints=constraints,priorprob="default")
BFreg2 <- BFregUpdate(BFreg1,lm1)








