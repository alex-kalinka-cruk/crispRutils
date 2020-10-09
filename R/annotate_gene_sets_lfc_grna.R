#' annotate_gene_sets_lfc_grna
#' 
#' Annotates gene sets in log2FC gRNA data and combines control and treatment data for plotting.
#' 
#' @param qc An object of class `fgcQC` (`fgcQC::QC_fgc_crispr_data`).
#' @return An object of class `lfc_gene_sets` which contains a data frame with a `type` column referring to control/treatment and a `gene_set` column.
#' @export
#' @importFrom dplyr mutate case_when
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
annotate_gene_sets_lfc_grna <- function(qc){
  if(!inherits(qc,"fgcQC")) stop(paste("expecting an object of class 'fgcQC', got",class(qc)))
  
  data <- qc$log2FC$control_vs_plasmid.gRNA %>%
    dplyr::mutate(type = "Control") %>%
    rbind(qc$log2FC$treatment_vs_plasmid.gRNA %>%
            dplyr::mutate(type = "Treatment")) %>%
    dplyr::mutate(gene_set = dplyr::case_when(gene %in% crispr_gene_sets$essential$pan_cancer_Sanger ~ "pan_cancer_Sanger",
                                              gene %in% crispr_gene_sets$essential$hart_essential ~ "hart_essential",
                                              gene %in% crispr_gene_sets$essential$moderately_negative ~ "moderately_negative",
                                              gene %in% crispr_gene_sets$essential$weakly_negative ~ "weakly_negative",
                                              gene %in% crispr_gene_sets$essential$hart_nonessential ~ "hart_nonessential",
                                              TRUE ~ "other"),
                  gene_set = factor(gene_set, levels = c("pan_cancer_Sanger","hart_essential","moderately_negative",
                                                         "weakly_negative","hart_nonessential"))) %>%
    dplyr::filter(gene_set != "other")
  ret <- list(data = data, type = "gRNA")
  class(ret) <- "lfc_gene_sets"
  return(ret)
}
