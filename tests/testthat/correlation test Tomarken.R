

getwd()
schiz <- read.table("/Users/jorismulder/downloads/Schiz_stuff/sz.txt",header=T)

dyn.load("/Users/jorismulder/surfdrive/paper BFpack/Rcode/BCT/test BCT/bct_prior_BFpack.dll")
dyn.load("/Users/jorismulder/surfdrive/paper BFpack/Rcode/BCT/test BCT/bct_continuous_final_BFpack.dll")

devtools::install_github("jomulder/BFpack")
library(BFpack)
schiz$Group <- as.factor(schiz$Group)
lm1 <- lm(cbind(Im_fr,De_fr,WM_let,Cat_fl,Fas,Rat) ~ -1 + Group,schiz)
BFmlm(lm1,parameter="correlation")

hypothesis = "
  Im_fr_with_De_fr_in_Group1 > Im_fr_with_De_fr_in_Group2 &
  Im_fr_with_De_fr_in_Group1 > Im_fr_with_WM_let_in_Group2 &
  Im_fr_with_De_fr_in_Group1 > Im_fr_with_Cat_fl_in_Group2 &
  Im_fr_with_De_fr_in_Group1 > Im_fr_with_Fas_in_Group2 &
  Im_fr_with_De_fr_in_Group1 > Im_fr_with_Rat_in_Group2 &

  De_fr_with_WM_let_in_Group1 > De_fr_with_WM_let_in_Group2 &
  De_fr_with_Cat_fl_in_Group1 > De_fr_with_Cat_fl_in_Group2 &
  De_fr_with_Fas_in_Group1 > De_fr_with_Fas_in_Group2 &
  De_fr_with_Rat_in_Group1 > De_fr_with_Rat_in_Group2 &

  WM_let_with_Cat_fl_in_Group1 > WM_let_with_Cat_fl_in_Group2 &
  WM_let_with_Fas_in_Group1 > WM_let_with_Fas_in_Group2 &
  WM_let_with_Rat_in_Group1 > WM_let_with_Rat_in_Group2 &

  Cat_fl_with_Fas_in_Group1 > Cat_fl_with_Fas_in_Group2 &
  Cat_fl_with_Rat_in_Group1 > Cat_fl_with_Rat_in_Group2 &

  Fas_with_Rat_in_Group1 > Fas_with_Rat_in_Group2;

Im_fr_with_De_fr_in_Group1 = Im_fr_with_De_fr_in_Group2 &
  Im_fr_with_De_fr_in_Group1 = Im_fr_with_WM_let_in_Group2 &
  Im_fr_with_De_fr_in_Group1 = Im_fr_with_Cat_fl_in_Group2 &
  Im_fr_with_De_fr_in_Group1 = Im_fr_with_Fas_in_Group2 &
  Im_fr_with_De_fr_in_Group1 = Im_fr_with_Rat_in_Group2 &

  De_fr_with_WM_let_in_Group1 = De_fr_with_WM_let_in_Group2 &
  De_fr_with_Cat_fl_in_Group1 = De_fr_with_Cat_fl_in_Group2 &
  De_fr_with_Fas_in_Group1 = De_fr_with_Fas_in_Group2 &
  De_fr_with_Rat_in_Group1 = De_fr_with_Rat_in_Group2 &

  WM_let_with_Cat_fl_in_Group1 = WM_let_with_Cat_fl_in_Group2 &
  WM_let_with_Fas_in_Group1 = WM_let_with_Fas_in_Group2 &
  WM_let_with_Rat_in_Group1 = WM_let_with_Rat_in_Group2 &

  Cat_fl_with_Fas_in_Group1 = Cat_fl_with_Fas_in_Group2 &
  Cat_fl_with_Rat_in_Group1 = Cat_fl_with_Rat_in_Group2 &

  Fas_with_Rat_in_Group1 = Fas_with_Rat_in_Group2

  "

library(nlme)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
lm1 <- lmList(cbind(Im_fr,De_fr) ~ 1 | Group,schiz)


