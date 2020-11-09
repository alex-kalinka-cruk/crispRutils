# Helper functions.

# LogFC scatterplot for Treatment-vs-Plasmid vs Control-vs-Plasmid.
.plot_plasmid_fc <- function(x, multimappers, label){
  data <- x$mageck_ctrl_treat_vs_plasmid
  if(!is.null(multimappers)){
    # Remove multi-mappers and controls from the plot.
    data %<>%
      dplyr::filter(!id %in% multimappers & id != "Control")
  }
  paired <- data %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(data = data,
                        ggplot2::aes(neg.lfc.Control, neg.lfc.Treatment, color = `-log10_FDR`)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::ggtitle("Essential (vs Plasmid) logFC KO vs Parental")
  
  if("Gain-Ess" %in% x$mageck_treat_vs_ctrl_ess_annot$type){
    new_ess <- (x$mageck_treat_vs_ctrl_ess_annot %>%
                  dplyr::filter(type == "Gain-Ess"))$id
    if(label){
      paired <- paired +
        ggplot2::geom_text(data = data %>% 
                             dplyr::filter(id %in% new_ess),
                           ggplot2::aes(neg.lfc.Control, neg.lfc.Treatment, label = id, hjust=0, vjust=0))
    }else{
      paired <- paired +
        ggplot2::geom_point(data = data %>% 
                              dplyr::filter(id %in% new_ess),
                            ggplot2::aes(neg.lfc.Control, neg.lfc.Treatment, fill = "Gain-Ess"),
                            pch = 21, color = "NA")
    }
  }
  print(paired)
}

# Volcano plot for Treatment-vs-Control.
.plot_treatment_volcano <- function(x, multimappers, label){
  data <- x$mageck_treat_vs_ctrl_ess_annot
  if(!is.null(multimappers)){
    # Remove multi-mappers and controls from the plot.
    data %<>%
      dplyr::filter(!id %in% multimappers & id != "Control")
  }
  
  volc <- data %>%
    ggplot2::ggplot(ggplot2::aes(log2FC,`-log10_FDR`)) +
    ggplot2::geom_point() +
    ggtitle("logFC(Treatment/Control)")
  
  if("Gain-Ess" %in% data$type){
    volc <- volc +
      ggplot2::geom_point(data = data %>%
                   dplyr::filter(type == "Gain-Ess"),
                 ggplot2::aes(log2FC,`-log10_FDR`, color = "Gain-Ess"))
    if(label){
      volc <- volc +
        ggplot2::geom_text(data = data %>% 
                  dplyr::filter(type == "Gain-Ess"),
                ggplot2::aes(log2FC, `-log10_FDR`, label = id, hjust = 0, vjust = 0))
    }
  }
  
  if("Loss-Ess" %in% data$type){
    volc <- volc + 
      ggplot2::geom_point(data = data %>%
                            dplyr::filter(type == "Loss-Ess"),
                          ggplot2::aes(log2FC,`-log10_FDR`, color = "Loss-Ess"))
  }
  
  print(volc)
}


#' plot.synthetic_lethal
#' 
#' Plot of Treatment-vs-plasmid against Control-vs-plasmid (logFC scatterplots), and volcano plot of Treatment-vs-Control. New essential genes and lost essential genes are highlighted, but only new essential genes are annotated with gene names.
#' 
#' @param x An object of class `synthetic_lethal`.
#' @param type Either `plasmid` or `treat`.
#' @param remove_multimappers Logical, whether to remove genes that have at least 1 guide that is a multi-mapper. Issues a warning if any of the genes appear in the new or lost essential groups.
#' @param label_hits Logical, whether to add gene names to plot, or simply colour them. Defaults to `TRUE`.
#' @param library_name Name of the library, e.g. "yusa_v3_human" - must match naming conventions used in the AZ-CRUK reference data repo.
#' @param library_type One of "n" (knock-out), "a" (activation), or "i" (interference).
#' @param library_annotation_version An integer giving the annotation version of the library.
#' @param ... Other arguments to pass to `plot`.
#' @return Renders a plot in the console.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_text aes geom_hline geom_vline ggtitle
#' @importFrom dplyr case_when filter
#' @importFrom magrittr %<>%
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
plot.synthetic_lethal <- function(x, type, remove_multimappers, 
                                  label_hits = TRUE,
                                  library_name = NULL, library_type = NULL,
                                  library_annotation_version = NULL, ...){
  if(!inherits(x,"synthetic_lethal")) stop(paste("expecting object of class 'synthetic_lethal', got:",class(x)))
  
  # Determine the multi-mapping genes.
  if(remove_multimappers){
    if(any(sapply(c(library_name, library_type, library_annotation_version),is.null)))
      stop("when 'remove_multimappers' is TRUE then library_name, library_type, and library_annotation_version must all be supplied")
    lt <- dplyr::case_when(library_type == "n" ~ "knockout",
                           library_type == "a" ~ "activation",
                           library_type == "i" ~ "interference",
                           TRUE ~ "unknown")
    if(lt == "unknown") stop(paste("library_type must be one of 'n', 'a', or 'i', got:",library_type))
    tryCatch(
      multimappers <- crispr_libs_annot_data[[lt]][[library_name]][[paste("v",library_annotation_version,sep="")]][["multimapping_guides"]][["genes"]],
      error = function(e) stop(paste("unable to extract multi-mapping genes for",library_name,library_type,library_annotation_version))
    )
    
    # Are any of the multi-mappers appearing in the new and lost essentials?
    new_ess_M <- intersect((x$mageck_treat_vs_ctrl_ess_annot %>%
      dplyr::filter(type == "Gain-Ess"))$id, multimappers)
    lost_ess_M <- intersect((x$mageck_treat_vs_ctrl_ess_annot %>%
                   dplyr::filter(type == "Loss-Ess"))$id, multimappers)
    
    if(length(new_ess_M) > 0)
      cat(paste("the following genes with multi-mapping gRNAs are in the 'Gain-Ess' group:\n",paste(new_ess_M, collapse=", ")))
    if(length(lost_ess_M) > 0)
      cat(paste("the following genes with multi-mapping gRNAs are in the 'Loss-Ess' group:\n",paste(lost_ess_M, collapse=", ")))
  }else{
    multimappers <- NULL
  }
  
  switch(type,
         plasmid = .plot_plasmid_fc(x, multimappers, label_hits),
         treat = .plot_treatment_volcano(x, multimappers, label_hits)
         )
}
