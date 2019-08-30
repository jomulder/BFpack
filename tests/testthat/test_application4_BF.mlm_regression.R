fmri.lm <- lm(cbind(Superficial,Middle,Deep) ~ Face + Vehicle, data=fmri)
constraints.fmri <- "Face_on_Deep = Face_on_Superficial = Face_on_Middle < 0 <
    Vehicle_on_Deep = Vehicle_on_Superficial = Vehicle_on_Middle;
    Face_on_Deep < Face_on_Superficial = Face_on_Middle < 0 < Vehicle_on_Deep =
    Vehicle_on_Superficial = Vehicle_on_Middle"
set.seed(123)
BF_fmri <- BF(fmri.lm, hypothesis = constraints.fmri)
#check results
matrixmeasures <- matrix(c(0.2662937, 0.1963028, 0.7248133, 0.9998789,  2.721856,  5.093553,  13.86392, 0.024443699,
                           0.3062633, 0.0251314, 5.1518462, 0.8251533, 16.821626, 32.833554, 552.31376, 0.973793185,
                           1.0000000, 1.0000000, 1.0000000, 1.0000000,  1.000000,  1.000000,   1.00000, 0.001763116),byrow=T,nrow=3)
expect_equivalent(
  BF_fmri$BFtable_confirmatory,matrixmeasures
)

