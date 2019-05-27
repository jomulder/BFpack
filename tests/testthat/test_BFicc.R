
# test BF testing on icc's / between-groups variances
library(bain)
library(BFpack)
# parse_hypothesis <- getFromNamespace("parse_hypothesis", "bain")
# make_RrList <- getFromNamespace("make_RrList", "BFpack")
# make_RrList2 <- getFromNamespace("make_RrList2", "BFpack")

# Use lmer-object as input
# remove(list=ls())
# library(lme4)
sleepstudy$category1 <- c(rep(1,30),rep(0,110),rep(1,40))
sleepstudy$category2 <- c(rep(0,30),rep(1,110),rep(0,40))
sleepstudy <- sleepstudy[sample(1:180,180,replace=F),]
# this randomly changes cluster indices without the respective outcome. Just for an example...
levels(sleepstudy$Subject) <- levels(sleepstudy$Subject)[sample(1:18,18,replace=F)]
fit1 <- lmer(Reaction ~ -1 + Days + category1 + category2 +
            (0 + category2 | Subject) + (0 + category1 | Subject), sleepstudy)

# only exploratory testing
BFtest1 <- BF(fit1)
# both exploratory and confirmatory testing
constraints="icc1=icc2>icc3>0;icc1>0 & icc2>icc3;icc1=0"
BFtest2 <- BF(fit1,hypothesis=constraints)
