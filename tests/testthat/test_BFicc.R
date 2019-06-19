#
# # test BF on icc's / between-groups variances on lmerMod object
#
# # library(bain)
# # library(BFpack)
# # # parse_hypothesis <- getFromNamespace("parse_hypothesis", "bain")
# # # make_RrList <- getFromNamespace("make_RrList", "BFpack")
# # # make_RrList2 <- getFromNamespace("make_RrList2", "BFpack")
# #
# # # Use lmer-object as input
# # # remove(list=ls())
# # # library(lme4)
# # sleepstudy$category1 <- c(rep(1,30),rep(0,110),rep(1,40))
# # sleepstudy$category2 <- c(rep(0,30),rep(1,110),rep(0,40))
# # sleepstudy <- sleepstudy[sample(1:180,180,replace=F),]
# # # this randomly changes cluster indices without the respective outcome. Just for an example...
# # levels(sleepstudy$Subject) <- levels(sleepstudy$Subject)[sample(1:18,18,replace=F)]
# # fit1 <- lmer(Reaction ~ -1 + Days + category1 + category2 +
# #             (0 + category2 | Subject) + (0 + category1 | Subject), sleepstudy)
# # fit1 <- lmer(Reaction ~ -1 + Days + category1 + category2 +
# #                (1 | Subject), sleepstudy)
# # fit1 <- lmer(Reaction ~ 1 + Days +
# #                (0 + category2 | Subject) + (0 + category1 | Subject), sleepstudy) # fixed intercepts
# # fit1 <- lmer(Reaction ~ -1 + Days + category1 + category2 +
# #                (Days | Subject), sleepstudy)
# # # only exploratory testing
# # BFtest1 <- BF(fit1)
# # # both exploratory and confirmatory testing
# # constraints="icc_category2=icc_category1=0;icc_category2>icc_category1>0"
# # BFtest2 <- BF(fit1,hypothesis=constraints)
# #
# #
# #
# #
# # ############# TIMMS Analyse JP Fox
# # load("/Users/jorismulder//Dropbox/JP-Joris/paper testing intraclass correlation/Real Data/TIMMS2015/data.Rdata")
# #
# #
# # XF <- cbind(as.numeric(sdata3m$ITSEX)-1,sdata3m$TOTWGT,sdata3m$year*(as.numeric(sdata3m$ITSEX)-1),sdata3m$lln)
# # yX <- cbind(sdata3m$math,sdata3m$group11a,sdata3m$group11b,sdata3m$group2a,sdata3m$group2b,sdata3m$group7a,
# #             sdata3m$group7b,sdata3m$group4a,sdata3m$group4b,XF) ## Complete Data Set
# # p <- 15
# #
# # ## 15 observations per school, 8 countries, ngroups = number of schools per country
# # ## Construct grouping variables for schools in each country c
# #
# # groupschool <- list()
# # for(ii in 1:8){
# #   groupschool[[ii]] <- as.factor(unlist(lapply(1:ngroups[ii],function(jj,p){rep(jj,p)},p=p)))
# # }
# #
# # ## Make Dataframe
# # yXTot <- data.frame(yX)
# #
# # names(yXTot) <- c("math","groupNL11","groupNL15","groupHR11","groupHR15","groupDE11","groupDE15","groupDK11","groupDK15",
# #                   "gender","weight","year","lln")
# # yXTot$groupschool <- unlist(groupschool)
# library(lme4)
# outlme <- lmer(math ~ -1 + gender + weight + yeargender + lln +
#                  groupNL11 + (0+groupNL11 | schoolID) +
#                  groupHR11 + (0+groupHR11 | schoolID) +
#                  groupDE11 + (0+groupDE11 | schoolID),data=timssICC)
# summary(outlme)
# BF(outlme)
#
#
# #
# #
# # ########
# #
# # load("C:\\Users\\JP-Su\\Dropbox\\JP-Joris\\paper testing intraclass correlation\\Real Data\\BF-lmer\\data-BFlmer.RData")
# # load("/Users/jorismulder/Dropbox/JP-Joris/data-BFlmer.RData")
# #
# # library(lme4)
# # ## for one year (2011) and countries (4 countries)
# timssICC_subset <- subset(timssICC, timssICC$groupNL11==1 | timssICC$groupHR11==1) #| timssICC$groupDE11==1 | timssICC$groupDK11==1)
# outlme1 <- lmer(math ~ -1 + gender + weight + lln
#                   + groupNL11 + (0+groupNL11 | schoolID)
#                   + groupHR11 + (0+groupHR11 | schoolID)
#                   # + groupDE11 + (0+groupDE11 | schoolID)
#                   # + groupDK11 + (0+groupDK11 | schoolID)
#                 ,data=timssICC_subset)
# summary(outlme1)
# BF(outlme1)
# BFout <- BF(outlme1,hypothesis="0<icc_groupNL11<icc_groupHR11;
#             0=icc_groupNL11=icc_groupHR11")
# # BFout <- BF(outlme1,hypothesis="icc_groupNL11<icc_groupHR11;
# #             icc_groupNL11=icc_groupHR11")
# #
# # BFout
# # BFout$estimates
#
#
#
#
#
#
