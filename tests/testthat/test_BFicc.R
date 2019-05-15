
# test BF testing on icc's / between-groups variances
library(bain)
library(BFpack)
parse_hypothesis <- getFromNamespace("parse_hypothesis", "bain")
make_RrList <- getFromNamespace("make_RrList", "BFpack")
make_RrList2 <- getFromNamespace("make_RrList2", "BFpack")

# Use lmer-object as input
sleepstudylist <- fm <- list()
sleepstudylist[[1]] <- sleepstudy[1:70,]
sleepstudylist[[2]] <- sleepstudy[71:120,]
sleepstudylist[[3]] <- sleepstudy[121:180,]
for(ca in 1:length(sleepstudylist)){
  # the fitted lmer-objects should always have fixed intercepts as well,
  # and only a random intercept as a random effect.
  fm[[ca]] <- lmer(Reaction ~ 1 + Days + (1 | Subject), sleepstudylist[[ca]])
}

# only exploratory testing
BFtest1 <- BFlmer(fm)
# both exploratory and confirmatory testing
constraints="icc1=icc2>icc3>0;icc1>0 & icc2>icc3;icc1=0"
BFtest2 <- BFlmer(fm,hypothesis=constraints)
