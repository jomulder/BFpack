
# test for correlation test on hetcor object
my_data <- mtcars[, c(1,3,4,5,6,7)]
res <- polycor::hetcor(my_data)
BF1 <- BF(res)
summary(BF1)
BF1 <- BF(res,hypothesis="disp_with_mpg < wt_with_mpg < 0;disp_with_mpg = wt_with_mpg < 0")
summary(BF1)
BF(res,hypothesis="disp_with_mpg < wt_with_mpg < 0;disp_with_mpg = wt_with_mpg < 0",prior=c(1,2,3))




