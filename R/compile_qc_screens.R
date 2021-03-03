#' compile_qc_screens
#' 
#' Compiles QC data from multiple screens.
#' 
#' @param path A path to a set of AZ-CRUK CRISPR screen results containing QC data (`QC_fgc.rds`).
#' @return A data frame.
#' @export
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
compile_qc_screens <- function(path){
  if(!dir.exists(path)) stop(paste("unable to find",path))
  qc_files <- list.files(path, recursive = T, pattern = "QC_fgc.rds$", full.names = T)
  if(length(qc_files) == 0) stop(paste("no 'QC_fgc.rds' files found in",path))
  
  ret <- NULL
  for(qc_file in qc_files){
    ret <- rbind(ret, readRDS(file=qc_file)$qc_metrics)
  }
  return(ret)
}
