

# cor1 <- cor_test_continuous(mtcars[,4:6],iter = 100)
# print(cor1$covmF)
# cat(cor1$covmF)
#
# cor1 <- cor_test(mtcars[,4:6],iter = 100)
# print(cor1$covmF)
# cat(cor1$covmF)

print("test bct")

set.seed(123)

Ntot <- 20
data <- mvtnorm::rmvnorm(Ntot,mean=c(0,0,0),sigma=1+diag(3))
K <- 1
P <- 3
numG <- 1
numcorr <- 3
BHat <- array(as.double(apply(data[,1:3],2,mean)),dim=c(numG,K,P))
sdHat <- array(as.double(apply(data[,1:3],2,sd)),dim=c(numG,P))
SigmaHat <- t(data[,1:3] - rep(1,Ntot)%*%t(BHat[1,1,]))%*%(data[,1:3] - rep(1,Ntot)%*%t(BHat[1,1,]))/(Ntot-1)
CHat <- array(as.double(diag(1/sdHat[1,]) %*% SigmaHat %*% diag(1/sdHat[1,])),dim=c(numG,P,P))
XtXi <- array(as.double(1/Ntot),dim=c(numG,K,K))
burnin <- 1e2
samsize0 <- 1e1
Ygroups <- array(0,dim=c(numG,Ntot,P))
Ygroups[1,,] <- round(data[,1:3],3)
Xgroups <- array(1,dim=c(numG,Ntot,K))
Njs <- as.double(Ntot)
sdMH <- array(as.double(apply(Ygroups[1,,],2,sd)/sqrt(nrow(Ygroups[1,,]))),dim=c(numG,P))
ordi <- Cat <- matrix(0,nrow=numG,ncol=P)
gLiuSab <- array(0,dim=c(samsize0,numG,P))
Njs <- matrix(as.double(Ntot),nrow=numG,ncol=1)

res <- .Fortran("estimate_bct_ordinal_test",
                P=as.integer(P),
                numG=as.integer(numG),
                Ntot=as.integer(Ntot),
                Ygroups=Ygroups,
                Cnugget=array(0,dim=c(P,P))
)

print(res$Cnugget)
print(res$Ygroups)
print(res$Ntot)


set.seed(123)
P <- 3
samsize <- 5000
Fisher <- 1
testm <- matrix(0,ncol=.5*P*(P-1),nrow=samsize)
res <-.Fortran("draw_ju",P = as.integer(P),
               drawscorr=testm,
               samsize=as.integer(samsize),
               numcorrgroup=as.integer(.5*P*(P-1)),
               Fisher=as.integer(Fisher),
               seed=as.integer( sample.int(1e6,1) )
)
print(head(res$drawscorr))
print(tail(res$drawscorr))


print("test cor_test_cont")
set.seed(123)
fit <- cor_test_continuous(BFpack::memory[,1:3])
print(head(fit$corrdraws[[1]][,2,1]))
print(tail(fit$corrdraws[[1]][,2,1]))


