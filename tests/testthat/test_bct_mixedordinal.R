

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
data1 <- mvtnorm::rmvnorm(Ntot,mean=c(0,0,0),sigma=1+diag(3))
qwe <- cor_test_test(as.data.frame(data1))
print(qwe)
print(qwe$corrdraws[[1]][100,1,])



print("prior draws")
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







