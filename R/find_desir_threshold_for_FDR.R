# Helper function.
.prep_4_mle <- function(data, desir_thresh){
  data %<>%
    dplyr::mutate(dplyr::across(-gene, function(x)  as.integer(x >= desir_thresh))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total = sum(dplyr::c_across(-gene))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(total) & gene != "Control")
  return(data)
}


#' find_desir_threshold_for_FDR
#' 
#' Finds the Desirability threshold that corresponds to a given FDR % for **shared** hits.
#' 
#' @param data A data frame of Desirability scores produced by `crispRutils::join_desirability_scores`.
#' @param start_desir A numeric value giving the desirability value to start at.
#' @param stop_desir A numeric value giving the desirability value to finish at.
#' @param fdr_thresh A numeric value giving the FDR threshold to aim for. Defaults to 0.1 (10% FDR).
#' @param step_size A numeric value giving the step size for sliding from the start to finish Desirability thresholds. Defaults to 0.05.
#' 
#' @return A data frame of FDR estimates at each Desirability threshold.
#' @importFrom dplyr %>% mutate across rowwise ungroup filter
#' @importFrom magrittr %<>%
#' @importFrom perept EM.perept
#' @export
#' @md
#' @author Alex T. Kalinka \email{alex.kalinka@@cancer.org.uk}
find_desir_threshold_for_FDR <- function(data, start_desir, stop_desir, fdr_thresh = 0.1, step_size = 0.05){
  tryCatch({
    desir_thresh <- start_desir
    finish <- FALSE
    track_perf <- NULL
    while(!finish){
      cat("Desirability Threshold:",desir_thresh,"\n")
      td <- .prep_4_mle(data, desir_thresh)
      mle <- perept::EM.perept(td$total, ncol(data)-1)
      # Estimate of global FDR.
      num_TP <- mle$prevalence$theta_1*nrow(td)
      num_TN <- nrow(td)-num_TP
      num_FP <- mle$performance$p10*num_TN
      FDR <- num_FP/(num_FP+num_TP)
      # Estimate of FDR after excluding singletons.
      if(length(which(td$total>=2))==0){
        # No shared hits.
        FDR.comm <- 0
        prob_pos_2 <- 0
      }else{
        prob_two_false <- 1-mle$delta$delta_pos[which(td$total==2)[1]]
        prob_three_false <- ifelse(sum(td$total==3)>0, 1-mle$delta$delta_pos[which(td$total==3)[1]], 0)
        num_two_false <- round(length(which(td$total==2))*prob_two_false)
        num_three_false <- ifelse(sum(td$total==3)>0, round(length(which(td$total==3))*prob_three_false), 0)
        FDR.comm <- (num_two_false + num_three_false)/(length(which(td$total>=2)))
        prob_pos_2 <- mle$delta$delta_pos[which(td$total==2)[1]]
      }
      track_perf <- rbind(track_perf, data.frame(Desir_threshold = desir_thresh,
                                                 FPR = mle$performance$p10,
                                                 FDR.global = FDR,
                                                 FDR.common = FDR.comm,
                                                 Sensitivity = mle$performance$p11,
                                                 Prevalence = round(mle$prevalence$theta_1*nrow(td)),
                                                 Prob_pos.two_shared_hits = prob_pos_2,
                                                 Number_singletons = length(which(td$total==1)),
                                                 Number_common_hits = length(which(td$total>=2))))
      if(!is.null(stop_desir)){
        if(desir_thresh <= stop_desir){
          finish <- TRUE
          ret <- td %>%
            mutate(probability_hit = mle$delta$delta_pos)
        }
      }
      if(FDR.comm >= fdr_thresh){
        finish <- TRUE
        ret <- td %>%
          mutate(probability_hit = mle$delta$delta_pos)
      }else{
        desir_thresh <- desir_thresh - step_size
      }
    }
  },
  error = function(e) stop(paste("unable to find Desirability threshold for given FDR:",e))
  )
  return(track_perf)
}
