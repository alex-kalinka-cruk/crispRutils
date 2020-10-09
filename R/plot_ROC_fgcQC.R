#' plot_ROC_fgcQC
#'
#' Plot ROC curves for sets of essential genes.
#'
#' @param data A data frame containing the following columns: `True_Positive_Rate`, `False_Positive_Rate`, and `gene_set`.
#' @param title The plot title.
#' @param print_auroc Logical indicating whether a table of AUROC values should be printed to the console. Defaults to `TRUE`.
#'
#' @return A ROC plot for sets of genes.
#' @importFrom ggplot2 ggplot aes geom_line geom_abline ggtitle
#' @export
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
plot_ROC_fgcQC <- function(data, title, print_auroc = TRUE){
  tryCatch({
    pl <- ggplot2::ggplot(data, ggplot2::aes(False_Positive_Rate, True_Positive_Rate, color = gene_set)) +
      ggplot2::geom_line() +
      ggplot2::geom_abline(linetype = "dashed") +
      ggplot2::ggtitle(paste("ROC curves:",title))
    print(pl)
    # Return a table of AUROCs.
    if(print_auroc)
      print(crispRutils::make_AUROC_df(data))
  },
  error = function(e) stop(paste("unable to plot ROC curves:",e))
  )
}
