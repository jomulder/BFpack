

# testing coefficients in multivariate normal model
lm1 <- lm(cbind(mpg,cyl,hp) ~ disp + wt, data = mtcars)
BF1 <- BF(lm1)
summary(BF1)
# tests on same predictor on different DVs
BF1 <- BF(x=lm1,hypothesis="disp_on_mpg>disp_on_cyl>disp_on_hp>0;disp_on_mpg=disp_on_cyl=disp_on_hp=0")
summary(BF1)
# tests on different predictors on same DVs
BF1 <- BF(lm1,hypothesis="disp_on_mpg>wt_on_mpg;disp_on_mpg=wt_on_mpg;disp_on_mpg<wt_on_mpg")
summary(BF1)
# tests on different predictors on different DVs
BF(lm1,hypothesis="disp_on_mpg>disp_on_cyl>0;disp_on_mpg=wt_on_cyl=0")


# testing correlations in multivariate normal model
lm1 <- lm(cbind(mpg,cyl,hp) ~ disp + wt, data = mtcars)
BF(lm1,parameter="correlation")
# (dummy) hypotheses on the correlations
BF(lm1,parameter="correlation",hypothesis="cyl_with_mpg<hp_with_mpg<hp_with_cyl<0;
   cyl_with_mpg<hp_with_mpg<hp_with_cyl<0")
# (dummy) hypotheses on the correlations where a predictor variable is a factor that reflects group specification
mtcars$vsfac <- as.factor(mtcars$vs)
lm1 <- lm(cbind(mpg,cyl,hp) ~ -1 + disp + wt + vsfac, data = mtcars)
BF(lm1,parameter="correlation",hypothesis="cyl_with_mpg_in_vsfac0 < cyl_with_mpg_in_vsfac1 &
                                           hp_with_mpg_in_vsfac0 < hp_with_mpg_in_vsfac1 &
                                           hp_with_cyl_in_vsfac0 < hp_with_cyl_in_vsfac1")

