#' read_mageck_gene_results
#'
#' Reads in all Mageck gene-level results (`*gene_summary.txt`) for a single screen put through the AZ-CRUK CRISPR pipeline.
#'
#' @param path A path to AZ-CRUK pipeline output for a single screen.
#' @return An object of class `mageck_gene`.
#' @export
#' @importFrom rlang sym
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
read_mageck_gene_results <- function(path){
  if(!dir.exists(path)) stop(paste("unable to find",path))
  tryCatch({
    files <- list.files(path, recursive = T, full.names = T, pattern = "gene_summary.txt$")
    ret <- list(essential_genes = NA)
    for(file in files){
      name <- gsub("^.*?\\/mageck\\/(\\S+?)\\.gene_summary.txt$","\\1",file)
      data <- read.delim(file, stringsAsFactors = F)
      ret[[rlang::sym(name)]] <- data
    }
  },
  error = function(e) stop(paste("unable to read Mageck gene-level results:",e))
  )
  class(ret) <- "mageck_gene"
  return(ret)
}
