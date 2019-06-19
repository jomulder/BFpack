# # Test class coxph
# library(BFpack)
# library(survival)
#
# fit <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) +
#              cluster(id), bladder1)
# BF(fit)
# BF(fit, "rx > size = number")
