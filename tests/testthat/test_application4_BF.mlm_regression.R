fmri.lm <- lm(cbind(Superficial,Middle,Deep) ~ Face + Vehicle, data=fmri)
constraints.fmri <- "Face_on_Deep = Face_on_Superficial = Face_on_Middle < 0 <
    Vehicle_on_Deep = Vehicle_on_Superficial = Vehicle_on_Middle;
    Face_on_Deep < Face_on_Superficial = Face_on_Middle < 0 < Vehicle_on_Deep =
    Vehicle_on_Superficial = Vehicle_on_Middle"
set.seed(123)
BF_fmri <- BF(fmri.lm, hypothesis = constraints.fmri)
#check results
matrixmeasures <- matrix(c(0.2741948, 0.19630283, 0.7186756, 0.9998622,  2.621041,  5.093468,  13.35019, 0.023006434,
  0.2967301, 0.02525591, 5.1451176, 0.8243133, 17.339388, 32.638428, 565.93036, 0.975270261,
  1.0000000, 1.00000000, 1.0000000, 1.0000000,  1.000000,  1.000000, 1.00000, 0.001723304),byrow=T,nrow=3)

expect_equivalent(
  BF_fmri$BFtable_confirmatory,round(matrixmeasures,7)
)

