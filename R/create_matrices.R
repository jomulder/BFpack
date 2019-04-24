#This code was written by Caspar van Lissa.

parse_hypothesis <- function(varnames, hyp){
  names_est <- varnames
  # Clean varnames and hyp, to turn parameter names into legal object names
  # Might be better to do this elsewhere
  hyp <- rename_function(hyp)
  varnames <- rename_function(varnames)
  
  # Check if varnames occur in hyp.
  hyp_params <- params_in_hyp(hyp)
  
  # Find full varnames using partial matching
  match_names <- charmatch(hyp_params, varnames)
  # Any varnames that do not match?
  if(anyNA(match_names)){
    stop("Some of the parameters referred to in the 'hypothesis' do not correspond to parameter names of object 'x'.\n  The following parameter names in the 'hypothesis' did not match any parameters in 'x': ",
         paste(hyp_params[is.na(match_names)], collapse = ", "),
         "\n  The parameters in object 'x' are named: ",
         paste(varnames, collapse = ", "))
  }
  if(any(match_names == 0)){
    stop("Some of the parameters referred to in the 'hypothesis' matched multiple parameter names of object 'x'.\n  The following parameter names in the 'hypothesis' matched multiple parameters in 'x': ",
         paste(hyp_params[match_names == 0], collapse = ", "),
         "\n  The parameters in object 'x' are named: ",
         paste(varnames, collapse = ", "))
  }
  # Substitute partial names in hypothesis with full names
  for(par in 1:length(hyp_params)){
    hyp <- gsub(hyp_params[par], varnames[match_names[par]], hyp)
  }
  # Substitute parameter names with full names
  hyp_params <- varnames[match_names]
  
  # Start parsing hypothesis here
  
  legal_varnames <- sapply(hyp_params, grepl, pattern = "^[a-zA-Z\\.][a-zA-Z0-9\\._]{0,}$")
  if(!all(legal_varnames)){
    stop("Could not parse the names of the 'estimates' supplied to parse_hypothesis(). Estimate names must start with a letter or period (.), and can be a combination of letters, digits, period and underscore (_).\nThe estimates violating these rules were originally named: ", paste("'", names_est[!legal_varnames], "'", sep = "", collapse = ", "), ".\nAfter parsing, these parameters are named: ", paste("'", hyp_params[!legal_varnames], "'", sep = "", collapse = ", "), call. = FALSE)
  }
  
  # Currently disabled: Is it a problem if there are parameters that don't occur in the hypothesis?
  # if(FALSE & any(!varnames %in% hyp_params)){
  #   message(
  #     "Some of the parameters in your model are not referred to in the 'hypothesis'. Make sure that this is what you intended.\n  Your hypothesis is: ",
  #     hyp,
  #     "\n  The parameters that are not mentioned are: ",
  #     paste(varnames[which(!varnames %in% hyp_params)], collapse = ", ")
  #   )
  # }
  # Ultimately, it would be best to check the reverse as well. This would also
  # allow us to make a string like "(a|b|c)" with the varnames, and insert it
  # into the regular expressions below. That would be faster than the complex
  # expressions currently used, and less prone to breakage.
  # varnames_in_hyp <- sapply(varnames, grepl, x = hyp)
  # if (FALSE) {
  #
  # }
  # End check
  hyp <- gsub("\\s", "", hyp)
  if(grepl("[><=]{2,}", hyp)) stop("Do not use combined comparison signs e.g., '>=' or '=='")
  
  original_hypothesis <- hyp_list <- strsplit(hyp, ";")[[1]]
  hyp_list <- lapply(hyp_list, function(x){ strsplit(x, "&")[[1]]})
  hyp_list <- lapply(hyp_list, function(x){unlist(lapply(x, expand_compound_constraints))})
  hyp_list <- lapply(hyp_list, function(x){unlist(lapply(x, expand_parentheses))})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, flip_inequality)})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, constraint_to_equation)})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, order_terms)})
  
  n_constraints <- as.vector(sapply(hyp_list, function(x){c(sum(grepl("=", x)), sum(grepl("[<>]", x)))}))
  nec_vec <- as.vector(sapply(hyp_list, function(x){sum(grepl("=", x))}))
  hyp_mat <- lapply(1:length(hyp_list), function(i){
    if(n_constraints[((i-1)*2)+1] == 0){
      ERr <- NULL
    } else {
      ERr <- t(sapply(hyp_list[[i]][grep("=", hyp_list[[i]])],
                      constraint_to_row, varnames = varnames))
    }
    if(n_constraints[((i-1)*2)+2] == 0){
      IRr <- NULL
    } else {
      IRr <- t(sapply(hyp_list[[i]][grep("[<>]", hyp_list[[i]])],
                      constraint_to_row, varnames = varnames))
    }
    rbind(ERr, IRr)
  })
  
  list(hyp_mat = hyp_mat, n_constraints = n_constraints, n_ec = nec_vec, original_hypothesis = original_hypothesis)
}

#' Create (in)equality constraint matrices
#'
#' Parses a character string describing a set of informative hypotheses,
#' and returns (in)equality constraint matrices
#'
#' Informative hypotheses specified as a character string by "hyp" should
#' adhere to the following simple syntax: \itemize{ \item Competing hypotheses
#' are separated by ";".  Thus, "a=b;a>b" means that H1: a=b, and H2: a>b.
#' \item Each individual hypothesis consists of a (series of) (in)equality
#' constraint(s). Every single (in)equality constraint is of the form "R1*mu1 +
#' R2*mu2+... = r", where capital Rs refer to numeric scaling constants, must
#' refer to the names of parameters in the model, and the lower case r refers
#' to a constant. Standard mathematical simplification rules apply; thus,
#' "R1*mu1 = R2*mu2" is equivalent to "R1*mu1 - R2*mu2 = 0".  \item Multiple
#' unrelated constraints within one hypothesis can be chained by "&". Thus,
#' "a=b&c=d" means that H1: a=b AND c=d.  \item Multiple related constraints
#' within one hypothesis can be chained by repeating the (in)equality operators
#' "=", "<", or ">". Thus, "a<b<c" means that H1: a < b AND b < c.  \item
#' Parameters can be grouped by placing them in a parenthesized, comma
#' separated list. Thus, "(a,b)>c" means that H1: a > c AND b > c.  Similarly,
#' "(a,b)>(c,d)" means that H1: a > c AND b > c AND b > c AND b > d.  }
#'
#' @aliases create_matrices create_matrices
#' @param varnames A character (vector of characters), containing names of
#' variables used in the hypotheses.  %Object of class \code{\link{lm}}, from
#' which the model parameters are extracted.
#' @param hyp A character string, containing a hypothesis (see Details).
#' @return A pair of named matrices for every hypothesis specified in the
#' \code{hyp} argument; one matrix named ERr, specifying equality constraints,
#' and one matrix named IRr, specifying inequality constraints.
#' @author Caspar van Lissa
#' @keywords internal utilities
create_matrices <- function(varnames, hyp){
  if(is.null(varnames)) stop("Please input proper linear model object")
  hyp <- gsub("\\s", "", hyp)
  if(grepl("[><=]{2,}", hyp)) stop("Do not use combined comparison signs e.g., '>=' or '=='")
  
  hyp_list <- strsplit(hyp, ";")[[1]]
  hyp_list <- lapply(hyp_list, function(x){ strsplit(x, "&")[[1]]})
  hyp_list <- lapply(hyp_list, function(x){unlist(lapply(x, expand_compound_constraints))})
  hyp_list <- lapply(hyp_list, function(x){unlist(lapply(x, expand_parentheses))})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, flip_inequality)})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, constraint_to_equation)})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, order_terms)})
  
  hyp_list <- unlist(lapply(hyp_list, function(x){
    ERr <- x[grep("=", x)]
    IRr <- x[grep("[<>]", x)]
    if(length(ERr) == 0){
      ERr <- NULL
    } else {
      ERr <- t(sapply(ERr, constraint_to_row, varnames = varnames))
      colnames(ERr)[ncol(ERr)]<- "="
    }
    if(length(IRr) == 0){
      IRr <- NULL
    } else {
      IRr <- t(sapply(IRr, constraint_to_row, varnames = varnames))
      colnames(IRr)[ncol(IRr)]<- ">"
    }
    list(ERr = ERr, IRr = IRr)
  }), recursive = FALSE)
  
  names(hyp_list) <- paste0(names(hyp_list), rep(1:(length(hyp_list)/2), each = 2))
  hyp_list
}







#' Expand compound constraints
#'
#' Takes a compound constraint, with multiple (in)equality operators, and
#' expands it into simple constraints.
#'
#'
#' @param hyp Character. An (in)equality constraint
#' @return A character vector with one element for each simple constraint
#' @keywords internal
expand_compound_constraints <- function(hyp){
  equality_operators <- gregexpr("[=<>]", hyp)[[1]]
  if(length(equality_operators) > 1){
    string_positions <- c(0, equality_operators, nchar(hyp)+1)
    return(sapply(1:(length(string_positions)-2), function(pos){
      substring(hyp, (string_positions[pos]+1), (string_positions[pos+2]-1))
    }))
  } else {
    return(hyp)
  }
}





#' Expand parentheses
#'
#' Takes a constraint with parentheses containing a "vector" of parameter
#' labels, and recursively expands the parenthesized "vectors".
#'
#'
#' @param hyp Character. An (in)equality constraint
#' @return A character vector with one element for each simple constraint
#' @keywords internal
expand_parentheses <- function(hyp){
  parenth_locations <- gregexpr("[\\(\\)]", hyp)[[1]]
  if(!parenth_locations[1] == -1){
    if(length(parenth_locations) %% 2 > 0) stop("Not all opening parentheses are matched by a closing parenthesis, or vice versa.")
    expanded_contents <- strsplit(substring(hyp, (parenth_locations[1]+1), (parenth_locations[2]-1)), ",")[[1]]
    if(length(parenth_locations) == 2){
      return(paste0(substring(hyp, 1, (parenth_locations[1]-1)), expanded_contents, substring(hyp, (parenth_locations[2]+1), nchar(hyp))))
    } else {
      return(apply(
        expand.grid(expanded_contents, expand_parentheses(substring(hyp, (parenth_locations[2]+1), nchar(hyp)))),
        1, paste, collapse = ""))
    }
  } else {
    return(hyp)
  }
}





#' Flip inequality
#'
#' Takes a constraint and flips elements on both sides of any "<"
#' operators so that only the ">" operator is used to define inequality
#' constraints.
#'
#'
#' @param hyp Character. An (in)equality constraint
#' @return Character
#' @keywords internal
flip_inequality <- function(hyp){
  if(grepl("<", hyp)){
    loc <- gregexpr("<", hyp)[[1]][1]
    return(paste0(substring(hyp, (loc+1)), ">", substring(hyp, 1, (loc-1))))
  } else {
    return(hyp)
  }
}





#' Constraint to equation
#'
#' Formats a constraint as an equation that can be evaluated. Adds scalars
#' to all parameters in the constraint, and adds a XXXconstant parameter to any
#' constants. Also adds "*" and "+" operators where necessary for evaluation.
#'
#'
#' @param hyp Character. An (in)equality constraint
#' @return Character
#' @keywords internal
constraint_to_equation <- function(hyp){
  # If scalar comes after variable name, move it in front
  hyp <- gsub("([\\.a-zA-Z][a-zA-Z0-9_\\.]{0,})\\*(\\d+)", "\\2*\\1", hyp, perl = TRUE)
  
  # When the string starts with a word, OR when a word is not preceded by
  # a number or *-sign, replace the "word" with "1*word"
  # hyp <- gsub("(^|(?<![\\*\\d]))([a-zA-Z][a-zA-Z0-9_]{0,})", "1*\\2", hyp, perl = TRUE)
  # Gu: add if not being preceded by any of .a-zA-Z0-9 , then ignore. otherwise x.A becomes x.1*A
  hyp <- gsub("(^|(?<![\\*\\d\\.a-zA-Z0-9_]))([a-zA-Z\\.][a-zA-Z0-9_\\.]{0,})", "1*\\2", hyp, perl = TRUE)
  
  # If a number starts with a period, add a leading zero
  # hyp <- gsub("(^|(?<!\\d))(\\.[0-9]+)", "0\\2", hyp, perl = TRUE)
  # Gu: add if being preceded by any of a-zA-Z0-9_, then ignore. Otherwise x.1 becomes x0.1
  hyp <- gsub("(^|(?<![a-zA-Z0-9_]))(\\.[0-9]+)", "0\\2", hyp, perl = TRUE)
  
  # When a number is directly followed by a word, insert a * sign
  # hyp <- gsub("(\\d)(?=[a-zA-Z][a-zA-Z0-9_]{0,})", "\\1*", hyp, perl = TRUE)
  # Gu: but not precededed by any of a-zA-Z_. Otherwise, x1x becomes x1*x
  hyp <- gsub("(?<![a-zA-Z_\\.])(\\d)(?=[a-zA-Z][a-zA-Z0-9_]{0,})", "\\1*", hyp, perl = TRUE) #####
  
  # When the string starts with a floating point number, OR
  # when a floating point number is not preceded by + or -, add a +
  # hyp <- gsub("(^|(?<![+-]))([0-9]{1,}\\.[0-9]{0,})", "+\\2", hyp, perl = TRUE)
  # Gu: add a condition: followed by a *, because otherwise x0.1x becomes x+0.1x.
  hyp <- gsub("(^|(?<![+-]))([0-9]{1,}\\.[0-9]{0,}(?=\\*))", "+\\2", hyp, perl = TRUE)
  
  # When the string starts with an integer number, OR
  # when an integer number is not preceded by + or -, add a +
  # hyp <- gsub("(^|(?<![\\.0-9+-]))([0-9]{1,})(?!\\.)", "+\\2", hyp, perl = TRUE)
  # Gu: add a condition: followed by a *, because otherwise x1x becomes x+1x.
  hyp <- gsub("(^|(?<![\\.0-9+-]))([0-9]{1,})(?!\\.)(?=\\*)", "+\\2", hyp, perl = TRUE)
  
  # Gu: To complement the above two, for a single number (constant), add a +
  hyp <- gsub("(?<=[=<>]|^)(\\d+|(\\d+\\.\\d+))(?=[=+-<>]|$)", "+\\1", hyp, perl = TRUE)
  
  # When a number is followed by =, >, <, +, or -, or is the last character of
  # the string, add "*XXXconstant" to the number
  # hyp <- gsub("(\\d)((?=[=<>+-])|$)", "\\1*XXXconstant", hyp, perl = TRUE)
  # Gu: but not preceded by word, number, period and  underline
  gsub("(?<![a-zA-Z0-9_\\.])(\\d+|(\\d+\\.\\d+))((?=[=<>+-])|$)", "\\1*XXXconstant", hyp, perl = TRUE)
}






#' Order terms
#'
#' Moves all terms on the right hand side of a constraint, which has been
#' formatted as an equation, to the left hand side.
#'
#'
#' @param hyp Character. An (in)equality constraint formatted as equation.
#' @return Character
#' @keywords internal
order_terms <- function(hyp){
  eq_location <- gregexpr("[=<>]", hyp)[[1]]
  rhs <- substring(hyp, eq_location+1, nchar(hyp))
  rhs <- gsub("(?=[+-])", "-", rhs, perl = TRUE)
  paste0(rhs, substring(hyp, 1, eq_location))
}





#' Constraint to row
#'
#' Evaluate a constraint, which has been formatted as an equation, with
#' all terms moved to the left hand side, and return a single row of a
#' (in)equality constraint matrix.
#'
#'
#' @param hyp Character. An (in)equality constraint formatted as equation,
#' with all terms on the left hand side.
#' @return Numeric vector.
#' @keywords internal
constraint_to_row <- function(varnames, hyp){
  e <- new.env()
  objects <- c(varnames, "XXXconstant")
  invisible(sapply(objects, assign, value = 0, envir = e))
  constraint_expression <- parse(text = substring(hyp, 1, nchar(hyp)-1))
  e[["XXXconstant"]] <- -1
  equal_to <- eval(constraint_expression, envir = e)
  e[["XXXconstant"]] <- 0
  
  c(sapply(varnames, function(x){
    e[[x]] <- 1
    constant <- eval(constraint_expression, envir = e)
    e[[x]] <- 0
    constant
  }), equal_to)
}

params_in_hyp <- function(hyp){
  params_in_hyp <- trimws(unique(strsplit(hyp, split = "[ =<>,\\(\\);&\\*+-]+", perl = TRUE)[[1]]))
  params_in_hyp <- params_in_hyp[!sapply(params_in_hyp, grepl, pattern = "^[0-9]*\\.?[0-9]+$")]
  params_in_hyp[grepl("^[a-zA-Z]", params_in_hyp)]
}



rename_function <- function(text){
  fulltext <- paste(text, collapse = "")
  new_names <- names_est <- text
  #if(grepl("[\\(\\)]", fulltext)){
  #  text <- gsub("\\(", "___O___", text)
  #  text <- gsub("\\)", "___C___", text)
  #}
  text[text == "(Intercept)"] <- "Intercept"
  if(grepl(":", fulltext)){
    text <- gsub(":", "___X___", text)
  }
  
  if(grepl("mean of ", fulltext)){
    text <- gsub("mean of the differences", "difference", text)
    text <- gsub("mean of ", "", text)
  }
  
  # If any variables are subsetted from data.frames: remode the df part of the name
  remove_df <- sapply(text, grepl, pattern = "[\\]\\$]+", perl = TRUE)
  if(any(remove_df)){
    text[remove_df] <- sapply(text[remove_df], function(x){
      tmp_split <- strsplit(x, "[\\]\\$]+", perl = TRUE)[[1]]
      if(length(tmp_split)==1){
        x
      } else {
        tail(tmp_split, 1)
      }
    })
  }
  
  text
}

