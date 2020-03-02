#' make_CRISPR_analysis_comparisons
#' 
#' Create a data frame with all combinations of desired CRISPR sample combinations.
#' 
#' @param samples_1 First set of samples (e.g. all baseline samples).
#' @param samples_2 Second set of samples (e.g. all control samples).
#' @param num_reps The number of replicates to take in each combination.
#' 
#' @return A data frame containing all combinations of 2 replicates from the first and second set of samples.
#' @import combinat combn
#' @import dplyr mutate select rowwise ungroup
make_CRISPR_analysis_comparisons <- function(samples_1, samples_2, num_reps = 2){
  c1 <- unlist(apply(combinat::combn(samples_1, num_reps),2,function(x) paste(x,collapse=".")))
  c2 <- unlist(apply(combinat::combn(samples_2, num_reps),2,function(x) paste(x,collapse=".")))
  ret <- data.frame(samples_2 = rep(c2,3)) %>%
    dplyr::mutate(samples_1 = rep(c1, each=3)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(comp1.1 = gsub("^(.*?)\\..*$","\\1",samples_1),
                  comp1.2 = gsub("^.*?\\.(.*)$","\\1",samples_1),
                  comp2.1 = gsub("^(.*?)\\..*$","\\1",samples_2),
                  comp2.2 = gsub("^.*?\\.(.*)$","\\1",samples_2)) %>%
    dplyr::select(-samples_1,-samples_2) %>%
    dplyr::ungroup()
  return(ret)
}
