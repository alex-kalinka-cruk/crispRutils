#' add_PRROC
#' 
#' Adds precision and recall columns to a data frame containing TP and a score metric (e.g. Bayes factor).
#' 
#' @param data A data frame containing a 'TP' column.
#' @param score_col A character string naming a column containing a score metric.
#' 
#' @return A data frame.
#' @import dplyr mutate filter select
#' @export
add_PRROC <- function(data, score_col){
  if(!"TP" %in% colnames(data))
    stop("expecting to find a 'TP' column in 'data'")
  if(!score_col %in% colnames(data))
    stop(paste("expecting to find a",score_col,"column in 'data'"))
  sc <- sym(score_col)
  
}