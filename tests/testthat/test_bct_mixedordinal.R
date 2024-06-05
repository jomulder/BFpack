

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
samsize0 <- 1e3
Ygroups <- array(0,dim=c(numG,Ntot,P))
Ygroups[1,,] <- round(data[,1:3],3)
Xgroups <- array(1,dim=c(numG,Ntot,K))
Njs <- as.double(Ntot)
sdMH <- array(as.double(apply(Ygroups[1,,],2,sd)/sqrt(nrow(Ygroups[1,,]))),dim=c(numG,P))
ordi <- Cat <- matrix(0,nrow=numG,ncol=P)
gLiuSab <- array(0,dim=c(samsize0,numG,P))
Njs <- matrix(as.double(Ntot),nrow=numG,ncol=1)

res <- .Fortran("estimate_bct_ordinal_test",
                postZmean=matrix(0,nrow=numcorr,ncol=1),
                postZcov=matrix(0,nrow=numcorr,ncol=numcorr),
                P=as.integer(P),
                numcorr=as.integer(numcorr),
                K=as.integer(K),
                numG=as.integer(numG),
                BHat=round(BHat,3),
                sdHat=round(sdHat,3),
                CHat=round(CHat,3),
                XtXi=XtXi,
                samsize0=as.integer(samsize0),
                burnin=as.integer(burnin),
                Ntot=as.integer(Ntot),
                Njs=as.integer(Njs),
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
                gLiuSab=gLiuSab,
                seed=as.integer( 121 ),
                nuggetscale = round(.995,3)
)
print(res$BDrawsStore[1,1,,])
print(res$Xgroups)
print(res$postZmean)
print(res$CDrawsStore[1,1,,])
cat(res$CDrawsStore[1:2,1,,])
#plot(res$CDrawsStore[,1,1,3])
postZmean_test <- c(.451,0.187,.448)
test_that("BF.cor_test exploratory hypotheses on correlations correctly evaluated", {
  expect_equivalent(
    postZmean_test,round(c(res$postZmean),3), tolerance = .1
  )
})

res <- cor_test_continuous(as.data.frame(data))
print(res$corrdraws[[1]][1,1:3,1:3])
cat(res$corrdraws[[1]][1,1:3,1:3])


