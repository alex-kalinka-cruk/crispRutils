#' calculate_library_width
#'
#' Calculates the library width of a plasmid sample, as defined in Imkeller et al. (2020) - 90th percentile divided by the 10th percentile.
#'
#' @param counts A vector of integer counts.
#'
#' @return A numeric value.
#' @export
#' @importFrom stats quantile
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
#' @md
#' @references
#' * Imkeller, K. et al. 2020. gscreend: modelling asymmetric count ratios in CRISPR screens to decrease experiment size and improve phenotype detection. *Genome Biol* **21**, 53.
calculate_library_width <- function(counts){
  tryCatch({
    lwd <- stats::quantile(counts, 0.9)/stats::quantile(counts, 0.1)
    names(lwd) <- NULL
    return(lwd)
  },
  error = function(e) stop(paste("unable to calculate library width:",e))
  )
}
