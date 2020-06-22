#' make_AUROC_df
#'
#' Makes an AUROC data frame.
#'
#' @param data A data frame containing the following columns: `True_Positive_Rate`, `False_Positive_Rate`, and `gene_set`.
#'
#' @return A data frame of AUROCs.
#' @importFrom dplyr mutate rename arrange %>%
#' @export
make_AUROC_df <- function(data){
  tryCatch({
    data_auroc <- data %>%
      dplyr::mutate(TP = !GENE %in% crispr_gene_sets$essential$hart_nonessential) %>%
      fgcQC::calc_AUC(score_col = "BF", group = "gene_set") %>%
      dplyr::rename(AUROC = AUC) %>%
      dplyr::arrange(dplyr::desc(AUROC))
  },
  error = function(e) stop(paste("unable to make AUROC data frame:",e))
  )
  return(data_auroc)
}
