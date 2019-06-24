

timssICC_subset <- timssICC[(timssICC$groupNL11==1)+(timssICC$groupHR11==1)>0,]
outlme1 <- lme4::lmer(math ~ -1 + gender + weight + lln
                  + groupNL11 + (0+groupNL11 | schoolID)
                  + groupHR11 + (0+groupHR11 | schoolID)
                ,data=timssICC_subset)
summary(outlme1)
BF1 <- BF(outlme1,hypothesis="0<groupNL11<groupHR11")
summary(BF1)







