#' plot_PrRc_fgcQC
#'
#' Plot Precision-Recall curves for sets of essential genes.
#'
#' @param data A data frame containing the following columns: `Recall`, `Precision`, and `gene_set`.
#' @param title The plot title.
#' @param print_auprrc Logical indicating whether a table of AUROC values should be printed to the console. Defaults to `TRUE`.
#'
#' @return A ROC plot for sets of genes.
#' @importFrom ggplot2 ggplot aes geom_line geom_hline ggtitle
#' @export
plot_PrRc_fgcQC <- function(data, title, print_auprrc = TRUE){
  tryCatch({
    pl <- ggplot2::ggplot(data, ggplot2::aes(Recall, Precision, color = gene_set)) +
      ggplot2::geom_line() +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed") +
      ggplot2::ggtitle(paste("Precision-Recall curves:",title))
    print(pl)
    # Return a table of AUPrRc.
    if(print_auprrc)
      print(crispRutils::make_AUPrRc_df(data))
  },
  error = function(e) stop(paste("unable to plot Precision-Recall curves:",e))
  )
}
