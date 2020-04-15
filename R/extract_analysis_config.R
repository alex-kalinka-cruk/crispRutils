#' extract_config_comparisons
#'
#' Extract the sample comparisons from an analysis config as part of pre-pipeline validation checks.
#'
#' @param file Path to an analysis config `json`` file.
#' @param out Either `df`` or `html`. Defaults to `df`.
#'
#' @return A data frame containing the sample comparisons.
#' @import dplyr mutate select
#' @import formattable formattable formatter style
#' @import jsonlite fromJSON
#' @import magrittr %>%
extract_config_comparisons <- function(file, out = "df"){
  if(!out %in% c("df","html")) stop("'out' must be one of 'df' or 'html")
  if(!file.exists(file)) stop(paste("unable to find",file))

  tryCatch(jlist <- jsonlite::fromJSON(file),
           error = function(e) stop("unable to process file",file,":",e))
  samps <- jlist$samples %>%
    dplyr::select(slx_id,name,replicate,indexes,label)
  if(out == "df"){
    return(jlist)
  }else{
    formattable(samps, list(
      label = formatter("span",
                        style = x ~ style(color = "red"))
    ))
  }
}
