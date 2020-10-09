#' add_essential_genes
#' 
#' Retrieves lists of essential genes for all comparisons in objects of class `mageck_gene` and `bagel_gene`.
#' 
#' @param x An object of class `mageck_gene` or `bagel_gene`.
#' @param mageck_fdr A value between 0 and 1 giving the Mageck FDR cutoff for essentiality. Defaults to 0.1, for 10 percent FDR.
#' @param bagel_thresh A numeric value giving the Bagel Bayes Factor threshold. if `NULL` (default), it's taken from the `bagel_gene` object.
#' @return An object of class `mageck_gene` or `bagel_gene` depending on the class of input argument `x`.
#' @export
#' @importFrom dplyr %>% filter
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
add_essential_genes <- function(x, mageck_fdr = 0.1, bagel_thresh = NULL){
  if(!inherits(x,"mageck_gene") && !inherits(x,"bagel_gene"))
    stop(paste("expecting an object of class 'mageck_gene' or 'bagel_gene', got:",class(x)))
  if(is.null(bagel_thresh) && inherits(x,"bagel_gene")){
    bagel_thresh <- x$bagel_bayes_factor_essentiality_threshold
  }
  tryCatch({
    x$essential_genes <- list()
    if(inherits(x,"mageck_gene")){
      if("Control_vs_Plasmid" %in% names(x)){
        x$essential_genes$Control_vs_Plasmid <- (x$Control_vs_Plasmid %>%
          dplyr::filter(neg.fdr < mageck_fdr))$id
      }
      if("Treatment_vs_Plasmid" %in% names(x)){
        x$essential_genes$Treatment_vs_Plasmid <- (x$Treatment_vs_Plasmid %>%
                                                   dplyr::filter(neg.fdr < mageck_fdr))$id
      }
    }else{
      if("Control_vs_Plasmid" %in% names(x)){
        x$essential_genes$Control_vs_Plasmid <- (x$Control_vs_Plasmid %>%
                                                   dplyr::filter(BF > bagel_thresh))$GENE
      }
      if("Treatment_vs_Plasmid" %in% names(x)){
        x$essential_genes$Treatment_vs_Plasmid <- (x$Treatment_vs_Plasmid %>%
                                                   dplyr::filter(BF > bagel_thresh))$GENE
      }
    }
  },
  error = function(e) stop(paste("unable to retrieve essential genes from object of class",class(x),":",e))
  )
  return(x)
}
