

fit1 <- aov(yield ~ -1 + K + P:K, npk)
fit1 <- lm(yield ~ -1 + K + P:K, npk)
model.matrix(fit1)
summary(fit1)
fit1$xlevels
fit1$assign
fit1$effects
fit1$rank
fit1$call

x <- fit1
Xmat <- model.matrix(x)
Ymat <- model.matrix(x)%*%x$coefficients + x$residuals




fit1 <- aov(yield ~ N + K + P + block, npk)
summary(fit1)
model.matrix(fit1)
fit1 <- aov(yield ~ N + K + P + block - 1, npk)
summary(fit1)
model.matrix(fit1)
fit1 <- aov(yield ~ P + K + N + block, npk)
summary(fit1)
model.matrix(fit1)
fit1 <- lm(yield ~ P + K + N, npk)
summary(fit1)
fit1 <- lm(yield ~ -1 + P + K + N, npk)
summary(fit1)
t.test(yield ~ P,npk,equal.var=T)
t.test(yield ~ K,npk,equal.var=T)

attr(fit1$terms,"order")
attr(fit1$terms,"factors")
#check if intercept is present
attr(fit1$terms,"intercept")
#check class of variables
attr(fit1$terms,"dataClasses")


?aov
fit1 <- aov(yield ~ block + N, npk)
summary(fit1)
model.matrix(fit1)
fit1 <- aov(yield ~ block, npk)
summary(fit1)
fit1 <- aov(yield ~ block + N*P*K, npk)
summary(fit1)
model.matrix(fit1)
unique(model.matrix(fit1))
tXX <- t(model.matrix(fit1))%*%model.matrix(fit1)
eigen(tXX)$values
Xmat <- model.matrix(fit1)
eigen(t(Xmat)%*%Xmat)$values

data(sleep)
plot(extra ~ group, data = sleep)
## traditional ANOVA gives a p value of 0.00283
summary(aov(extra ~ group + Error(ID/group), data = sleep))

## Build design matrix
group.column <- rep(1/c(-sqrt(2),sqrt(2)),each=10)
subject.matrix <- model.matrix(~sleep$ID - 1,data=sleep$ID)
## Note that we include no constant column
X <- cbind(group.column, subject.matrix)
library(BayesFactor)
bf.full <- nWayAOV(y = sleep$extra, X = X, gMap = c(0,rep(1,10)), rscale=c(.5,1))
exp(bf.full[['bf']])
bf.full2 <- lmBF(extra ~ group + ID, data = sleep, whichRandom = "ID")
bf.full2

fit1 <- aov(extra ~ group + Error(ID/group), data=sleep)
fit1 <- lm(extra ~ group, data=sleep)
model.matrix(fit1)





delivery.df = data.frame(
  Service = c(rep("Carrier 1", 15), rep("Carrier 2", 15),
              rep("Carrier 3", 15)),
  Destination = c(rep(c("Office 1", "Office 2", "Office 3",
                        "Office 4", "Office 5"), 9)),
  Time = c(15.23, 14.32, 14.77, 15.12, 14.05,
           15.48, 14.13, 14.46, 15.62, 14.23, 15.19, 14.67, 14.48, 15.34, 14.22,
           16.66, 16.27, 16.35, 16.93, 15.05, 16.98, 16.43, 15.95, 16.73, 15.62,
           16.53, 16.26, 15.69, 16.97, 15.37, 17.12, 16.65, 15.73, 17.77, 15.52,
           16.15, 16.86, 15.18, 17.96, 15.26, 16.36, 16.44, 14.82, 17.62, 15.04)
)
ggplot(delivery.df, aes(Time, Destination, colour = Service)) + geom_point()
delivery.mod1 = aov(Time ~ Destination*Service, data = delivery.df)
summary(delivery.mod1)





# NEXT  EXAMPLE
attach(InsectSprays)
data(InsectSprays)
str(InsectSprays)
mean(count[spray=="A"])
tapply(count, spray, mean)
tapply(count, spray, var)
tapply(count, spray, length)
boxplot(count ~ spray)
Photoperiod<-ordered(spray,levels=c("F","B","C","D","E","A"))
tapply(count,Photoperiod,mean)
is.factor(spray)
aov.out = aov(count ~ spray, data=InsectSprays)
model.matrix(aov.out)
summary(aov.out)
TukeyHSD(aov.out)
summary.lm(aov.out)


####### YET ANOTHER EXAMPLE
Input =("
Diet    Country  Weight_change
        A       USA      0.120
        A       USA      0.125
        A       USA      0.112
        A       UK       0.052
        A       UK       0.055
        A       UK       0.044
        B       USA      0.096
        B       USA      0.100
        B       USA      0.089
        B       UK       0.025
        B       UK       0.029
        B       UK       0.019
        C       USA      0.149
        C       USA      0.150
        C       USA      0.142
        C       UK       0.077
        C       UK       0.080
        C       UK       0.066
        ")
Data = read.table(textConnection(Input),header=TRUE)
Data$cov1 <- rnorm(nrow(Data))
Data$Country = factor(Data$Country,
                      levels=unique(Data$Country))
library(psych)
headTail(Data)
str(Data)
summary(Data)
rm(Input)
interaction.plot(x.factor     = Data$Country,
                 trace.factor = Data$Diet,
                 response     = Data$Weight_change,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")
model = lm(Weight_change ~ Country + Diet + Country:Diet,
           data = Data)
Anova(model,type = "II")


fit1 = aov(Weight_change ~ -1 + Country + Diet, data = Data)
summary(fit1)
model.matrix(fit1)
fit1 = aov(Weight_change ~ Country + Diet, data = Data)
summary(fit1)
model.matrix(fit1)
fit1$xlevels
fit1 = aov(Weight_change ~ -1 + Country*Diet, data = Data)
summary(fit1)
model.matrix(fit1)
fit1 = aov(Weight_change ~ Country*Diet, data = Data)
summary(fit1)
model.matrix(fit1)
fit1 = lm(Weight_change ~ Diet*Country, data = Data)
summary(fit1)
model.matrix(fit1)

fit1$xlevels

fit1$terms
attr(fit1$terms,"order")
attr(fit1$terms,"factors")
#check if intercept is present
attr(fit1$terms,"intercept")
#check class of variables
attr(fit1$terms,"dataClasses")

BF(fit1)




