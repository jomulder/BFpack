#Test glm class object
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
BF1 <- BF(glm.D93)
BF1$estimates
BF(glm.D93, hypothesis = "treatment2 = outcome3>0")


