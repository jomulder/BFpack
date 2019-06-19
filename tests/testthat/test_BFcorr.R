# # Example (dummy) correlation test
# dt <- as.data.frame(scale(mtcars))
#
# # exploratory testing of correlations of one correlation matrix
# model1 <- lm(cbind(mpg,cyl,hp) ~ disp + wt, data = dt)
# BFcorr1 <- BFcorr(model1,prior=NULL,constraints="exploratory",priorprob="default")
# #get posterior probabilities for the hypotheses
# BFcorr1$PHP
# # confirmatory test
# constraints <- "cyl_with_mpg > hp_with_mpg > hp_with_cyl; cyl_with_mpg < hp_with_mpg = cyl_with_hp ; cyl_with_mpg = hp_with_mpg = cyl_with_hp"
# BFcorr1 <- BFcorr(model1,prior=NULL,constraints=constraints,priorprob="default")
# #get posterior probabilities for the hypotheses
# BFcorr1$PHP
#
# # confirmatory testing of correlations in two correlation matrix
# model2 <- lm(cbind(mpg,cyl,hp) ~ carb + gear, data = dt)
# model <- list(model1,model2)
#
# constraints <- "cyl_with_mpg_gr1 > hp_with_mpg_gr1 > hp_with_cyl_gr1; cyl_with_mpg_gr1 < hp_with_mpg_gr1 = cyl_with_hp_gr1 ; cyl_with_mpg_gr2 = hp_with_mpg_gr2 = cyl_with_hp_gr1"
# #constraints <- "cyl_with_mpg > hp_with_mpg > hp_with_cyl; cyl_with_mpg < hp_with_mpg = cyl_with_hp ; cyl_with_mpg = hp_with_mpg = cyl_with_hp"
# BFcorr1 <- BFcorr(model,prior=NULL,constraints=constraints,priorprob="default")
# #get posterior probabilities for the hypotheses
# BFcorr1$PHP
# BFcorr2 <- BFcorrUpdate(BFcorr1,model) #update with same information new data contained in 'model' (just for the exercise)
# #get posterior probabilities for the hypotheses
# BFcorr2$PHP
