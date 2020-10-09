#' plot.lfc_gene_sets
#' 
#' Plot logFC boxplots for gene sets.
#' 
#' @param x An object of class `lfc_gene_sets` (`crispRutils::annotate_gene_sets_lfc_grna`).
#' @param notch Add notches to the boxplot. Defaults to `TRUE`.
#' @param ... Other arguments to pass to `plot`.
#' @return Plot.
#' @export
#' @importFrom ggplot2 ggplot geom_point geom_boxplot aes ggtitle facet_grid theme element_blank geom_hline
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
plot.lfc_gene_sets <- function(x, notch = TRUE, ...){
  if(!inherits(x,"lfc_gene_sets")) stop(paste("expecting an object of class 'lfc_gene_sets', got:",class(x)))
  
  pl <- x$data %>%
    ggplot2::ggplot(ggplot2::aes(gene_set, log2FC, color = gene_set)) +
    ggplot2::geom_point() +
    ggplot2::geom_boxplot(notch = notch) +
    ggplot2::facet_grid(~type) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank()) +
    ggplot2::ggtitle("log2FC by gene groups: gRNA level")
  print(pl)
}
