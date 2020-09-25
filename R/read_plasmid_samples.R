#' read_plasmid_samples
#'
#' Reads in multiple plasmid count samples and returns a single data frame.
#'
#' @param dir A character string giving a directory path containing one or more plasmid count files.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr %>% mutate select everything rename
#' @importFrom utils read.delim
read_plasmid_samples <- function(dir){
  if(!dir.exists(dir))
    stop(paste("unable to find",dir))
  tryCatch({
    files <- list.files(dir, recursive = F, full.names = T)
    ret <- NULL
    for(file in files){
      fn <- tail(unlist(strsplit(file,"/")),1)
      if(grepl("csv$",file)){
        name <- gsub("^(.*?)_plasmid_counts.csv$","\\1",fn)
        sname <- gsub("^(.*?)--.*$","\\1",name)
        td <- read.csv(file, stringsAsFactors = F) %>%
          dplyr::select(gRNA,Gene,Plasmid) %>%
          dplyr::mutate(library = name, batch = "Ong-et-al-2017", count_normalized = 2e7*Plasmid/sum(Plasmid)) %>%
          dplyr::rename(count = Plasmid)
        # Add GC content from library.
        if(any(grepl(sname,names(crispr_libs)))){
          clib <- crispr_libs[[names(crispr_libs)[grepl(sname,names(crispr_libs))]]]
          td %<>%
            dplyr::mutate(GC_percent = clib$GC_percent[match(gRNA,clib$CODE)])
        }else{
          td$GC_percent <- NA
        }
      }else if(grepl("txt$",file)){
        name <- gsub("^(.*?)_plasmid_counts.txt$","\\1",fn)
        lib <- gsub("^.*?-(.*)$","\\1",name)
        slib <- gsub("^(.*?)--.*$","\\1",lib)
        batch <- gsub("^(.*?)-.*$","\\1",name)
        td <- utils::read.delim(file, stringsAsFactors = F, header = F) %>%
          dplyr::mutate(library = lib, batch = batch, count_normalized = 2e7*V3/sum(V3)) %>%
          dplyr::rename(gRNA = V1, Gene = V2, count = V3)
        # Add GC content from library.
        if(any(grepl(slib,names(crispr_libs)))){
          clib <- crispr_libs[[names(crispr_libs)[grepl(slib,names(crispr_libs))]]]
          td %<>%
            dplyr::mutate(GC_percent = clib$GC_percent[match(gRNA,clib$CODE)])
        }else{
          td$GC_percent <- NA
        }
      }
      ret <- rbind(ret,td)
    }
  },
  error = function(e) stop(paste("unable to process plasmid count samples:",e))
  )
  return(ret %>%
           dplyr::select(library, batch, dplyr::everything()))
}
