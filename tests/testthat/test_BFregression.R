# test BF on lm object

lm1 <-  lm(wt ~ disp + drat + hp, mtcars)
BF1 <- BF(lm1)
BF(lm1,hypothesis="disp=drat=0;disp>drat>0")

#BF(x=lm1,hypothesis="disp<.2 & disp> -.2")

