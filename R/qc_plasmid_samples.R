#' qc_plasmid_samples
#'
#' Calculates several QC metrics for sequenced plasmid sample counts.
#'
#' @param data A data frame of plasmid counts, as produced by `crispRutils::read_plasmid_samples`.
#'
#' @return A data frame of summary QC metrics.
#' @export
#' @importFrom dplyr %>% summarise group_by ungroup n
#' @importFrom moments skewness
#' @importFrom energy dcor2d
qc_plasmid_samples <- function(data){
  tryCatch({
    qc <- data %>%
      dplyr::group_by(library, batch) %>%
      dplyr::summarise(Library_width = crispRutils::calculate_library_width(count),
                       Skewness = moments::skewness(count),
                       Sarles_bimodality = crispRutils::calculate_sarles_bimodality(count),
                       percent_zero_count_guides = 100*sum(count == 0)/dplyr::n(),
                       percent_less_30_count_guides = 100*sum(count < 30)/dplyr::n(),
                       distcorr_GC_content_counts = ifelse(!all(is.na(GC_percent)),
                                                           energy::dcor2d(count[!is.na(GC_percent)],
                                                                          GC_percent[!is.na(GC_percent)], type="U"),
                                                           NA)) %>%
      dplyr::ungroup()
  },
  error = function(e) stop(paste("unable to QC plasmid samples:",e))
  )
  return(qc)
}
