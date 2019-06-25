
#test for variances
vtest1 <- var_test(InsectSprays$count, InsectSprays$spray)
hypothesis <- "A<B<C;C<D>E;A=C=E"
BF(vtest,hypothesis)





