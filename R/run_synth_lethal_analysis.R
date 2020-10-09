#' run_synth_lethal_analysis
#' 
#' Performs a synthetic-lethal analysis for a single screen in which a genetically modified (GM; usually a gene knock-out) isogenic cell line is compared to the unmodified parental line: the GM line is considered a 'treatment' and the parental line is the 'control' arm (as in a drug-interaction screen).
#' 
#' @param path A path to AZ-CRUK pipeline output for a single 'treatment-control' screen.
#' @return An object of class `synthetic_lethal`.
#' @export
#' @importFrom dplyr %>% mutate case_when inner_join rowwise ungroup
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
run_synth_lethal_analysis <- function(path){
  if(!dir.exists(path)) stop(paste("unable to find",path))
  
  # 1. Read in Mageck and Bagel gene-level results and extract essential genes.
  mageck <- crispRutils::read_mageck_gene_results(path) %>%
    crispRutils::add_essential_genes()
  bagel <- crispRutils::read_bagel_gene_results(path) %>%
    crispRutils::add_essential_genes()
  
  # 2. Find the union of Mageck and Bagel essential genes.
  ess_union <- crispRutils::union_essentials(mageck, bagel)
  
  # 3. Define new essentials and lost essentials relative to 'treatment' cell line.
  new_ess <- setdiff(ess_union$Treatment_vs_Plasmid, ess_union$Control_vs_Plasmid)
  lost_ess <- setdiff(ess_union$Control_vs_Plasmid, ess_union$Treatment_vs_Plasmid)
  
  # 4. Annotate treat-ctrl mageck output with new and lost ess genes.
  tryCatch({
    sleth_mageck <- mageck$Treatment_vs_Plasmid %>%
      dplyr::mutate(type = dplyr::case_when(
        (neg.fdr < 0.1 & id %in% new_ess) ~ "Gain-Ess",
        (pos.fdr < 0.1 & id %in% lost_ess) ~ "Loss-Ess",
        TRUE ~ "Unchanged-Ess"),
        `-log10_FDR` = -log10(min(c(neg.fdr, pos.fdr))))
    },
    error = function(e) stop(paste("unable to annotate mageck output with gain and loss of essentiality:",e))
  )
  
  # 5. Join ctrl-plasmid and treat-plasmid for plotting purposes.
  tryCatch({
    ctrl_treat_mageck <- mageck$Control_vs_Plasmid %>%
      dplyr::inner_join(mageck$Treatment_vs_Plasmid, by = "id", suffix = c(".Control",".Treatment")) %>%
      dplyr::rowwise() %>%
      # For plotting we use the best p-val for each gene's essentiality.
      dplyr::mutate(`-log10_FDR` = -log10(min(c(neg.fdr.Control, neg.fdr.Treatment)))) %>%
      dplyr::ungroup()
  },
  error = function(e) stop(paste("unable to join ctrl-plasmid and treat-plasmid mageck output:",e))
  )
  
  ret <- list(mageck_res = mageck, bagel_res = bagel, essential_union = ess_union,
              new_essential = new_ess, lost_essential = lost_ess,
              mageck_ctrl_treat_vs_plasmid = ctrl_treat_mageck,
              mageck_treat_vs_ctrl_ess_annot = sleth_mageck)
  class(ret) <- "synthetic_lethal"
  return(ret)
}
