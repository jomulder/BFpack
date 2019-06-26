#Test glm class object
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
BF1 <- BF(glm.D93)
BF1$estimates
BF(x=glm.D93, hypothesis = "treatment2 > 2 & treatment2 < 3")
vec1 <- 1:2
names(vec1) <- c("t1","t2")
bain(vec1,"t1 < 3 & t1 > 2",Sigma=diag(2),n=20)

bain(x=c(mean1),hypothesis=constraints1,Sigma=Sigma1,n=999)




sesamesim$site <- as.factor(sesamesim$site)
# execute an analysis of variance using lm() which, due to the -1, returns
# estimates of the means per group
anov <- lm(postnumb~site-1,sesamesim)
# take a look at the estimated means and their names
coef(anov)
# set a seed value
set.seed(100)
# use the names to formulate and test hypotheses with bain
results <- bain(anov, "site1=0")
results <- bain(anov, "site1>2 & site1<3")



