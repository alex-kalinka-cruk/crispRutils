#' read_bagel_gene_results
#' 
#' Reads in all Bagel gene-level results (`*.bf`) for a single screen put through the AZ-CRUK CRISPR pipeline.
#' 
#' @param path A path to AZ-CRUK pipeline output for a single screen.
#' @param ess_thresh One of `orig`, `fgcQC` - indicating which essential gene sets are used for determining the Bagel threshold.
#' @return An object of class `bagel_gene`.
#' @export
#' @importFrom rlang sym
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
read_bagel_gene_results <- function(path, ess_thresh){
  if(!dir.exists(path)) stop(paste("unable to find",path))
  if(!ess_thresh %in% c("orig","fgcQC")) 
    stop(paste("expecting 'ess_thresh' to be either 'orig' or 'fgcQC', got:",ess_thresh))
  tfile_patt <- ifelse(ess_thresh == "orig","EssentialityThresh.tsv","EssentialityThresh_fgcQc.tsv")
  ret <- list(essential_genes = NA)
  tryCatch({
    tfile <- list.files(path, recursive = T, full.names = T, pattern = tfile_patt)
    if(length(tfile) != 1)
      stop(paste("found",length(tfile),"essentiality threshold files, expected 1:",paste(tfile, collapse=", ")))
    ret$essentiality_threshold_type <- ess_thresh
    ret$bagel_bayes_factor_essentiality_threshold <- read.delim(tfile, stringsAsFactors = F)$value[1]
    files <- list.files(path, recursive = T, full.names = T, pattern = "bf$")
    for(file in files){
      name <- gsub("^.*?\\/bagel\\/(\\S+?)\\.bf$","\\1",file)
      data <- read.delim(file, stringsAsFactors = F)
      ret[[rlang::sym(name)]] <- data
    }
  },
  error = function(e) stop(paste("unable to read Bagel gene-level results:",e))
  )
  class(ret) <- "bagel_gene"
  return(ret)
}
