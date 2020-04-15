#' make_CRISPR_comparison_table
#'
#' Create a data frame with all (two-replicate) combinations of CRISPR samples.
#'
#' @param samples_1 First set of samples (e.g. all plasmid samples).
#' @param samples_2 Second set of samples we will sub-sample from (e.g. all control samples).
#' @param num_reps The number of replicates to take in each combination.
#'
#' @return A data frame containing all combinations of 2 replicates from the first and second set of samples.
#' @import combinat combn
#' @import dplyr mutate select rowwise ungroup
make_CRISPR_comparison_table <- function(samples_1, samples_2, num_reps = 2){
  num_pl <- length(samples_1)
  c2 <- unlist(apply(combinat::combn(samples_2, num_reps),2,function(x) paste(x,collapse=".")))
  ret <- data.frame(samples_2 = c2) %>%
    dplyr::mutate(samples_1 = paste(samples_1,collapse=".")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(comp2.1 = gsub("^(.*?)\\..*$","\\1",samples_2),
                  comp2.2 = gsub("^.*?\\.(.*)$","\\1",samples_2)) %>%
    dplyr::select(samples_1,comp2.1,comp2.2) %>%
    dplyr::ungroup()
  return(ret)
}
