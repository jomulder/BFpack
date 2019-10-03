#' @importFrom utils getFromNamespace
#' @importFrom bain t_test bain
#' @importFrom stats approxfun coef complete.cases cov dbeta density dnorm dt lm median model.matrix
#' @importFrom stats nobs pchisq pnorm pt quantile rWishart rbeta rgamma rnorm rt sd setNames var vcov
parse_hypothesis <- getFromNamespace("parse_hypothesis", "bain")
constraint_to_equation <- getFromNamespace("constraint_to_equation", "bain")
constraint_to_row <- getFromNamespace("constraint_to_row", "bain")
expand_compound_constraints <- getFromNamespace("expand_compound_constraints", "bain")
expand_parentheses <- getFromNamespace("expand_parentheses", "bain")
flip_inequality <- getFromNamespace("flip_inequality", "bain")
order_terms <- getFromNamespace("order_terms", "bain")
params_in_hyp <- getFromNamespace("params_in_hyp", "bain")

