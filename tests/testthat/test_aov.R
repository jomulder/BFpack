data(tvprices)
aov1 <- aov(price ~ anchor*motivation,data=tvprices)
#names(aov1$coefficients)
BF1 <- BF(aov1,hypothesis="anchorrounded=motivationlow;
   anchorrounded<motivationlow;
   anchorrounded>motivationlow")



