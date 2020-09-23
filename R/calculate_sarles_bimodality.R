#' calculate_sarles_bimodality
#'
#' Calculates Sarle's bimodality metric, which provides a measure of a skewness when the given distribution is unimodal.
#'
#' @param data A vector of integer counts.
#'
#' @return A numeric value.
#' @export
#' @importFrom moments skewness kurtosis
calculate_sarles_bimodality <- function(data){
  tryCatch({
    return((moments::skewness(data)^2+1)/moments::kurtosis(data))
  },
  error = function(e) stop(paste("unable to calculate Sarle's bimodality metric:",e))
  )
}
