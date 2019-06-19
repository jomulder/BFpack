# #Test glm class object
# library(BFpack)
#
# counts <- c(18,17,15,20,10,20,25,13,12)
# outcome <- gl(3,1,9)
# treatment <- gl(3,3)
#
# glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
#
# BF(glm.D93)
# BF.glm(glm.D93, hypothesis = "treatment2 = treatment3")
