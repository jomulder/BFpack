#' @importFrom utils getFromNamespace
#' @importFrom bain t_test
parse_hypothesis <- getFromNamespace("parse_hypothesis", "bain")
constraint_to_equation <- getFromNamespace("constraint_to_equation", "bain")
constraint_to_row <- getFromNamespace("constraint_to_row", "bain")
expand_compound_constraints <- getFromNamespace("expand_compound_constraints", "bain")
expand_parentheses <- getFromNamespace("expand_parentheses", "bain")
flip_inequality <- getFromNamespace("flip_inequality", "bain")
order_terms <- getFromNamespace("order_terms", "bain")
params_in_hyp <- getFromNamespace("params_in_hyp", "bain")
rankMatrix <- getFromNamespace("rankMatrix","Matrix")
pmvt <- getFromNamespace("pmvt","mvtnorm")
dmvt <- getFromNamespace("dmvt","mvtnorm")
rref <- getFromNamespace("rref","pracma")
ginv <- getFromNamespace("ginv","MASS")

