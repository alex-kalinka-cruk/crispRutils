#' union_essentials
#' 
#' Returns the union of essentials called by Mageck and Bagel (from objects of class `mageck_gene` and `bagel_gene`).
#' 
#' @param mageck An object of class `mageck_gene`.
#' @param bagel An object of class `bagel_gene`.
#' @return An object of class `essential_union`.
#' @export
union_essentials <- function(mageck, bagel){
  if(!inherits(mageck,"mageck_gene")) stop(paste("expecting an object of class 'mageck_gene', got:",class(mageck)))
  if(!inherits(bagel,"bagel_gene")) stop(paste("expecting an object of class 'bagel_gene', got:",class(bagel)))
  if(is.na(mageck$essential_genes) || is.na(bagel$essential_genes))
    stop("please first run 'add_essential_genes()'")
  
  ret <- list()
  tryCatch({
    if("Control_vs_Plasmid" %in% names(mageck) & "Control_vs_Plasmid" %in% names(bagel)){
      ret$Control_vs_Plasmid <- union(mageck$essential_genes$Control_vs_Plasmid,
                                      bagel$essential_genes$Control_vs_Plasmid)
    }
    if("Treatment_vs_Plasmid" %in% names(mageck) & "Treatment_vs_Plasmid" %in% names(bagel)){
      ret$Treatment_vs_Plasmid <- union(mageck$essential_genes$Treatment_vs_Plasmid,
                                        bagel$essential_genes$Treatment_vs_Plasmid)
    }
  },
  error = function(e) stop(paste("unable to retrieve union of essentials:",e))
  )
  class(ret) <- "essential_union"
  return(ret)
}
