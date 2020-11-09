#' read_desirability_scores
#' 
#' Read in Desirability scores from the AZ-CRUK CRISPR analysis pipeline.
#' 
#' @param path A path to AZ-CRUK pipeline output for a single screen.
#' 
#' @return A data frame.
#' @export
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
read_desirability_scores <- function(path){
  if(!dir.exists(path)) stop(paste("unable to find",path))
  tryCatch({
    file <- list.files(path, recursive = T, full.names = T, pattern = "DesirabilityGeneRanking.tsv$")
    if(length(file) != 1)
      stop(paste("expecting one 'DesirabilityGeneRanking.tsv' file, got:",file))
    desir <- read.delim(file, stringsAsFactors = F)
  },
  error = function(e) stop(paste("unable to read Desirability scores:",e))
  )
  return(desir)
}
