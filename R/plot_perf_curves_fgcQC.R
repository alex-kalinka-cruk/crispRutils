#' plot_perf_curves_fgcQC
#'
#' Plots a ROC or Precision-Recall curve for sets of essential genes taking the `fgcQC` output R object as input.
#'
#' @param qc_fgcQC An object of class `fgcQC`.
#' @param type A character string naming the plot type. One of: `ROC` or `PrRc`.
#'
#' @return A ROC or Precision-Recall curve for sets of essential genes.
#' @export
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
plot_perf_curves_fgcQC <- function(qc_fgcQC, type){
  if(!inherits(qc_fgcQC,"fgcQC"))
    stop(paste("expecting input to be an object of class 'fgcQC', got an object of class:",class(qc_fgcQC)))

  tryCatch({
    if(type == "ROC"){
      for(i in 1:length(qc_fgcQC$bagel_ROC)){
        comp <- gsub("^bagel_(\\S+?)$","\\1",names(qc_fgcQC$bagel_ROC)[i])
        crispRutils::plot_ROC_fgcQC(qc_fgcQC$bagel_ROC[[i]], comp)
      }
    }else{
      for(i in 1:length(qc_fgcQC$bagel_PrRc)){
        comp <- gsub("^bagel_(\\S+?)$","\\1",names(qc_fgcQC$bagel_PrRc)[i])
        crispRutils::plot_PrRc_fgcQC(qc_fgcQC$bagel_PrRc[[i]], comp)
      }
    }
  },
  error = function(e) stop(paste("unable to plot performance curves:",e))
  )
}
