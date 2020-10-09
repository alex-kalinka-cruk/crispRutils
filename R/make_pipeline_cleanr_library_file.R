#' make_pipeline_cleanr_library_file
#'
#' Prepares a 'cleanr' library file (to be compatible with the AZ-CRUK CRISPR pipeline) from either an sgRNA sequence file in which genomic co-ordinates are encoded in sgRNA IDs, or a csv file format used for TechDev benchmarking libraries.
#'
#' @param file A character string naming a path to an sgRNA sequence file.
#' @param grna_id_column A character string naming the sgRNA column.
#' @param seq_column A character string naming the sequence column.
#' @param gene_column A character string naming the sequence column.
#' @param csv Logical indicating whether input is a `csv` file or not. Defaults to `FALSE`.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr %>% mutate select rowwise ungroup rename
#' @importFrom magrittr %<>%
#' @importFrom rlang sym !! :=
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
make_pipeline_cleanr_library_file <- function(file, grna_id_column, seq_column, gene_column, csv = FALSE){
  tryCatch({
    if(!csv){
      data <- read.delim(file, stringsAsFactors = F, strip.white = T)
      data %<>%
        dplyr::rowwise() %>%
        dplyr::mutate(!!rlang::sym(seq_column) := toupper(!!rlang::sym(seq_column)),
                      coords = gsub("^.*?_(ex.*)$","\\1",!!rlang::sym(grna_id_column)),
                      EXONE = gsub("^ex(\\d+?)_.*$","\\1",coords),
                      CHRM = gsub("^ex\\d+?_(\\d+?):.*$","\\1",coords),
                      STRAND = gsub("^ex\\d+?_\\d+?:\\S+?:(.).*$","\\1",coords),
                      STARTpos = gsub("^ex\\d+?_\\d+?:(\\d+?)-.*$","\\1",coords),
                      ENDpos = gsub("^ex\\d+?_\\d+?:\\d+?-(\\d+?):.*$","\\1",coords)) %>%
        dplyr::ungroup() %>%
        dplyr::rename(seq = !!rlang::sym(seq_column),
                      CODE = !!rlang::sym(grna_id_column),
                      GENES = !!rlang::sym(gene_column)) %>%
        dplyr::select(seq, CODE, GENES, EXONE, CHRM, STRAND, STARTpos, ENDpos)
    }else{
      data <- read.csv(file) %>%
        dplyr::mutate(!!rlang::sym(seq_column) := toupper(!!rlang::sym(seq_column)),
                      CHRM = "NA", STARTpos = "NA", ENDpos = "NA") %>%
        dplyr::select(sequence, !!rlang::sym(grna_id_column), !!rlang::sym(gene_column), target_exon,
                      CHRM, strand, STARTpos, ENDpos) %>%
        dplyr::rename(seq = sequence, CODE = !!rlang::sym(grna_id_column),
                      GENES = !!rlang::sym(gene_column),
                      EXONE = target_exon, STRAND = strand)
    }
  },
  error = function(e) stop(paste("unable to make cleanr library file:",e))
  )
  return(data)
}
