# Helper function.
.inner_join_desir <- function(data_list, type){
  nams <- names(data_list)
  # Negative or Positive.
  dcol <- rlang::sym(ifelse(type == "neg","TvCCvP_NegativeDesi","TvCCvP_PositiveDesi"))
  for(i in seq_along(data_list)){
    data_list[[i]] %<>%
      dplyr::select(gene, !!dcol) %>%
      dplyr::rename(!!rlang::sym(paste(as.character(dcol),".",nams[i],sep="")) := !!dcol)
  }
  ij <- data_list %>%
    purrr::reduce(dplyr::inner_join, by = "gene")
  return(ij)
}


#' join_desirability_scores
#' 
#' Performs inner join of desirability scores (from the AZ-CRUK CRISPR analysis pipeline) across 2 or more screens.
#' 
#' @param ... Two or more objects returned by `crispRutils::read_desirability_scores`.
#' @param names A character vector providing names for the desirability data frames in `...` - must be the same length as `...`.
#' @param type One of either `neg` or `pos`, indicating negative desirability scores or positive ones respectively.
#' 
#' @return An object of class `desirability`.
#' @importFrom rlang sym !! :=
#' @importFrom magrittr %<>%
#' @importFrom dplyr %>% select rename inner_join
#' @importFrom purrr reduce
#' @export
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
join_desirability_scores <- function(..., names, type){
  input <- list(...)
  if(length(input) < 2) stop("please provide 2 or more desirability score input data frames")
  if(length(names) != length(input)) stop("the number of names must be equal to the number of input data frames")
  if(!type %in% c("neg","pos")) stop(paste("expecting 'type' to be either 'neg' or 'pos', got:",type))
  names(input) <- names
  tryCatch({
    ret <- .inner_join_desir(input, type)
  },
  error = function(e) stop(paste("unable to join desirability score data frames:",e))
  )
  return(ret)
}
