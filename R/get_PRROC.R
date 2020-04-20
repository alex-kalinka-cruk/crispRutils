#' get_PRROC
#'
#' Returns precision and recall values for a data frame containing TP and a score metric (e.g. Bayes factor), or the AUPrRc.
#'
#' @param data A data frame containing a 'TP' column.
#' @param score_col A character string naming a column containing a score metric.
#' @param return A character string indicating what should be returned: if `curve` then the PR curve is added to the data frame, if `auc` then the AUCPrRc value is returned.
#'
#' @return A data frame or a single numeric value.
#' @import dplyr mutate filter select rename
#' @import PRROC pr.curve
#' @export
get_PRROC <- function(data, score_col, return = "curve"){
  if(!"TP" %in% colnames(data))
    stop("expecting to find a 'TP' column in 'data'")
  if(!score_col %in% colnames(data))
    stop(paste("expecting to find a",score_col,"column in 'data'"))
  sc <- sym(score_col)
  pos_scores <- c((data %>%
                   dplyr::filter(TP) %>%
                   dplyr::select(!!sc))[,1])
  neg_scores <- c((data %>%
                     dplyr::filter(!TP) %>%
                     dplyr::select(!!sc))[,1])
  pr <- PRROC::pr.curve(pos_scores, neg_scores, curve = T)
  if(return == "curve"){
    prc <- as.data.frame(pr$curve) %>%
      dplyr::rename(Recall = V1, Precision = V2)
    return(prc)
  }else{
    return(pr$auc.davis.goadrich)
  }
}
