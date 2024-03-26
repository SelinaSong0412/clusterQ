#' Identify interaction terms 
#' 
#' @description A function that identifying all interaction covariates with a target covariate in a model formula
#'
#' @param formula A formula that specifying the regression model 
#' @param target_covariate A string that indicating the name of the target covariate 
#' 
#' @return A vector of covariate names that interact with the target covariate
#' 
#' @export
#' 
#' @examples
#' model_formula <- formula(Y ~ A1 * X2 + A1:X3 + X4)
#' target_covariate <- "A1"
#' interacting_variables <- interact.terms(model_formula, target_covariate)

interact.terms <- function(formula, target_covariate) {
  
  formula_str <- deparse(formula, width.cutoff = 200)
  interacting_variables <- character()
  # Extract terms from the formula
  terms <- strsplit(formula_str, "~|\\s*(?![^(]*\\))\\s*[\\+\\-]\\s*", perl = TRUE)[[1]]
  # terms <- strsplit(formula_str, "\\+|\\-|\\~")
  for (term in terms) {
    if (grepl(target_covariate, term)) {
      if (grepl(":", term) | grepl("\\*", term)) {
        interaction_parts <- unlist(strsplit(term, split = "\\*|:", perl = TRUE))
        interaction_parts <- gsub("\\s+", "", interaction_parts) # Remove all spaces
        vars <- interaction_parts[interaction_parts != target_covariate]
        if (length(vars) > 0) {
          interacting_variables <- c(interacting_variables, vars)
        }
      }
    }
  }
  interacting_variables <- gsub("\\s+", "", interacting_variables) # Remove all spaces
  unique(interacting_variables)
}

