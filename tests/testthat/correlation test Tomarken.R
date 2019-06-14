

getwd()
schiz <- read.table("/Users/jorismulder/downloads/Schiz_stuff/sz.txt",header=T)

dyn.load("/Users/jorismulder/surfdrive/paper BFpack/Rcode/BCT/test BCT/bct_prior_BFpack.dll")
dyn.load("/Users/jorismulder/surfdrive/paper BFpack/Rcode/BCT/test BCT/bct_continuous_final_BFpack.dll")

devtools::install_github("jomulder/BFpack")
library(BFpack)
lm1 <- lm(cbind(Im_fr,De_fr,WM_let,Cat_fl,Fas,Rat) ~ 1,schiz)
BF(lm1,parameter="correlation")



