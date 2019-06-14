#test class glmerMod
library(BFpack)
library(lme4)

gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                  data = cbpp, family = binomial)

BF.glmerMod(gm1)
