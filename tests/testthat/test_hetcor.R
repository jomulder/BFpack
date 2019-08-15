data(fmri)
# test for correlation test on hetcor object
res <- hetcor(fmri[,3:5])
BF1 <- BF(res)
summary(BF1)
BF1 <- BF(res,hypothesis="(Middle_with_Superficial,Deep_with_Superficial,Deep_with_Middle) > 0;
          Middle_with_Superficial=Deep_with_Superficial=Deep_with_Middle= 0")
summary(BF1)




