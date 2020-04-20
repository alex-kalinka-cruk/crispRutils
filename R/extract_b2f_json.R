#' extract_b2f_json
#' 
#' Extracts flowcell-level QC metrics from the 'Stats.json' output from 'bcl2fastq2'.
#' 
#' @param path A character string giving a path to the 'Stats.json' output from 'bcl2fastq2'.
#' 
#' @import dplyr mutate select filter unnest
#' @import jsonlite fromJSON
#' @return A data frame
#' 
#' @export
extract_b2f_json <- function(path){
  ret <- list()
  all <- jsonlite::fromJSON(path) %>%
    dplyr::as_tibble()
  # Sanity check input.
  # Flowcell info.
  flowcell <- all[,1:3]
  # RunInfo (read config).
  run_info <- as.data.frame(all$ReadInfosForLanes$ReadInfos) %>%
    dplyr::mutate(RunId = flowcell$RunId)
  # Flowcell summary stats.
  summ <- all$ConversionResults %>%
    dplyr::select(-DemuxResults, -Undetermined) %>%
    dplyr::mutate(RunId = flowcell$RunId)
  # Sample summary stats.
  samps <- all$ConversionResults$DemuxResults[[1]] %>%
    dplyr::mutate(RunId = flowcell$RunId) %>%
    tidyr::unnest(c(IndexMetrics, ReadMetrics), 
                  names_repair="universal") %>%
    dplyr::mutate(MismatchCounts = MismatchCounts[,1],
                  Index_OneBaseMismatch_percent = 100*(NumberReads-MismatchCounts)/NumberReads,
                  Q30_bases_percent = 100*YieldQ30/Yield...8,
                  Average_base_quality = QualityScoreSum/Yield...8,
                  Trimmed_bases_percent = 100*TrimmedBases/Yield...8,
                  Sample_Representation = 100*NumberReads/summ$TotalClustersPF)
  not_demux_count <- all$ConversionResults$Undetermined$NumberReads
  summ %<>%
    dplyr::mutate(ReadsPF_percent = 100*TotalClustersPF/TotalClustersRaw,
                  Non_Demultiplexed_Reads_percent = 100*not_demux_count/summ$TotalClustersPF,
                  Q30_bases_samples_percent = 100*sum(samps$YieldQ30)/sum(samps$Yield...8))
  ret$flowcell <- flowcell
  ret$run_info <- run_info
  ret$summary <- summ
  ret$samples <- samps
  return(ret)
}
