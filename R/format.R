#' Get the variable labels for a vector of column names
#'
#' @param vars_names A data frame with columns "variable" (column name) and "label"
#' @param varnames A character vector containing column names.
#' @return The corresponding labels of the variables in \code{varnames}.
#'
#' @keywords Internal
#'
#'
get_label <- function(vars_names, varnames){
  out = data.frame(variable = varnames) %>%
    left_join(vars_names, by = "variable") %>%
    # If no label is found, just use the inputted name
    mutate(label = ifelse(is.na(label), variable, label))

  return(out$label)
}
