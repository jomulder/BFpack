set.seed(123)
#set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion")
lm6 <- lm(cbind(Im,Del,Wmn,Cat,Fas,Rat) ~ -1 + Group, data=memory)
BF6_cor <- BF(lm6,parameter="correlation", hypothesis=
    "Del_with_Im_in_GroupHC > Del_with_Im_in_GroupSZ &
    Del_with_Wmn_in_GroupHC > Del_with_Wmn_in_GroupSZ &
    Del_with_Cat_in_GroupHC > Del_with_Cat_in_GroupSZ &
    Del_with_Fas_in_GroupHC > Del_with_Fas_in_GroupSZ &
    Del_with_Rat_in_GroupHC > Del_with_Rat_in_GroupSZ &
    Im_with_Wmn_in_GroupHC > Im_with_Wmn_in_GroupSZ &
    Im_with_Cat_in_GroupHC > Im_with_Cat_in_GroupSZ &
    Im_with_Fas_in_GroupHC > Im_with_Fas_in_GroupSZ &
    Im_with_Rat_in_GroupHC > Im_with_Rat_in_GroupSZ &
    Wmn_with_Cat_in_GroupHC > Wmn_with_Cat_in_GroupSZ &
    Wmn_with_Fas_in_GroupHC > Wmn_with_Fas_in_GroupSZ &
    Wmn_with_Rat_in_GroupHC > Wmn_with_Rat_in_GroupSZ &
    Cat_with_Fas_in_GroupHC > Cat_with_Fas_in_GroupSZ &
    Cat_with_Rat_in_GroupHC > Cat_with_Rat_in_GroupSZ &
    Fas_with_Rat_in_GroupHC > Fas_with_Rat_in_GroupSZ")
#check results
expect_equivalent(
  round(BF6_cor$BFmatrix_confirmatory[1,2],3),6201.679
)

#
# install.packages("ggplot2")
# ggplot(attention, aes(x = group, y = accuracy)) +
#   geom_density_ridges(aes(fill = group)) +
#   scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))

