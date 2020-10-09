#' make_AUPrRc_df
#'
#' Makes an AUPrRc data frame.
#'
#' @param data A data frame containing the following columns: `Recall`, `Precision`, and `gene_set`.
#'
#' @return A data frame of AUPrRcs.
#' @importFrom dplyr group_by summarise ungroup arrange desc
#' @export
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
make_AUPrRc_df <- function(data){
  tryCatch({
    data_auprrc <- data %>%
      dplyr::group_by(gene_set) %>%
      dplyr::summarise(AUPrRc = AUPrRc[1],
                       Sensitivity_FDR_10pct = Sensitivity_FDR_10pct[1],
                       Sensitivity_FDR_5pct = Sensitivity_FDR_5pct[1],
                       .groups = 'drop') %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(AUPrRc))
  },
  error = function(e) stop(paste("unable to make AUPrRc data frame:",e))
  )
  return(data_auprrc)
}
