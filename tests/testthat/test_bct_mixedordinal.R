

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
BHat <- array(0,dim=c(numG,K,P))
BHat[1,1,] <- apply(data[,1:3],2,mean)
sdHat <- array(1,dim=c(numG,P))
sdHat[1,] <- apply(data[,1:3],2,sd)
SigmaHat <- t(data[,1:3] - rep(1,Ntot)%*%t(BHat[1,1,]))%*%(data[,1:3] - rep(1,Ntot)%*%t(BHat[1,1,]))/(Ntot-1)
CHat <- array(1,dim=c(numG,P,P))
CHat[1,,] <- diag(1/sdHat[1,]) %*% SigmaHat %*% diag(1/sdHat[1,])
XtXi <- array(1/Ntot,dim=c(numG,K,K))
burnin <- 1e2
samsize0 <- 1e2
Ygroups <- array(0,dim=c(numG,Ntot,P))
Ygroups[1,,] <- data[,1:3]
Xgroups <- array(1,dim=c(numG,Ntot,K))
Njs <- Ntot
sdMH <- array(apply(Ygroups[1,,],2,sd)/sqrt(nrow(Ygroups[1,,])),dim=c(numG,P))
ordi <- Cat <- matrix(c(0,0,0),nrow=numG,ncol=P)
gLiuSab <- array(0,dim=c(samsize0,numG,P))
Njs <- matrix(Ntot,nrow=numG,ncol=1)

res <- .Fortran("estimate_bct_ordinal",
                postZmean=matrix(0,nrow=numcorr,ncol=1),
                postZcov=matrix(0,nrow=numcorr,ncol=numcorr),
                P=as.integer(P),
                numcorr=as.integer(numcorr),
                K=as.integer(K),
                numG=as.integer(numG),
                BHat=BHat,
                sdHat=sdHat,
                CHat=CHat,
                XtXi=XtXi,
                samsize0=as.integer(samsize0),
                burnin=as.integer(burnin),
                Ntot=as.integer(Ntot),
                Njs_in=Njs,
                Xgroups=Xgroups,
                Ygroups=Ygroups,
                C_quantiles=array(0,dim=c(numG,P,P,3)),
                sigma_quantiles=array(0,dim=c(numG,P,3)),
                B_quantiles=array(0,dim=c(numG,K,P,3)),
                BDrawsStore=array(0,dim=c(samsize0,numG,K,P)),
                sigmaDrawsStore=array(0,dim=c(samsize0,numG,P)),
                CDrawsStore=array(0,dim=c(samsize0,numG,P,P)),
                sdMH=sdMH,
                ordinal_in=ordi,
                Cat_in=Cat,
                maxCat=as.integer(max(Cat)),
                gLiuSab=array(0,dim=c(samsize0,numG,P)),
                seed=as.integer( 121 ),
                nuggetscale = .995)
print(res$CDrawsStore[1:2,1,,])
cat(res$CDrawsStore[1:2,1,,])

