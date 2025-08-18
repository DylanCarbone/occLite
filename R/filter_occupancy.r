#' Flag and filter years with very wide confidence intervals
#'
#' Filters each species' occupancy Index to remove years where the absolute
#' confidence-interval width exceeds a threshold.
#'
#' @param occti_outputs A named list of occti outputs, where each element has
#'   a component \code{$Index} containing columns \code{Year}, \code{psiA_L},
#'   and \code{psiA_U}.
#' @param ci_width_threshold Numeric scalar. Years with \code{psiA_U - psiA_L}
#'   greater than or equal to this value are removed. Default is \code{0.8}.
#'
#' @details
#' For each species in \code{occti_outputs}, the function computes the CI width
#' as \code{psiA_U - psiA_L}, drops rows with CI width \eqn{\ge} \code{ci_width_threshold},
#' and writes a short message listing the removed years (if any).
#'
#' @return The same \code{occti_outputs} list, updated in place so that each
#'   element's \code{$Index} only contains years with CI width < \code{ci_width_threshold}.
#'   No additional objects are returned.
#'
#' @examples
#' # occti_outputs <- list(
#' #   "Species A" = list(Index = data.frame(Year=2000:2005, psiA_L=c(.1,.2,.3,.05,.1,.2), psiA_U=c(.6,.7,1,.95,.5,.55))),
#' #   "Species B" = list(Index = data.frame(Year=2001:2004, psiA_L=c(.05,.1,.2,.2), psiA_U=c(.9,.95,.7,.6)))
#' # )
#' # occti_outputs <- filter_occupancy(occti_outputs, ci_width_threshold = 0.8)
filter_occupancy <- function(occti_outputs, ci_width_threshold = 0.8, z_score_threshold = 3, min_years_post_filter = 3) {

    species_to_remove = c()

    for (species in names(occti_outputs)){

        occ_data <- occti_outputs[[species]]$Index

        if (nrow(occ_data) == 1) {
        message(paste("species", species, "has only one year of data. Removing species"))
        species_to_remove <- c(species_to_remove, species)
        next
        }

        occ_data_filtered <- occ_data %>%
            mutate(ci_width = psiA_U - psiA_L,
            mean_psiA = mean(psiA),
            z_score = abs((psiA - mean_psiA)/sd(psiA)),
            ci_width_flag = ifelse(ci_width > ci_width_threshold, TRUE, FALSE),
            z_score_flag = ifelse(z_score > z_score_threshold, TRUE, FALSE)
            )
        
        # Provide message for species-years dropped
        if(sum(occ_data_filtered$ci_width_flag) > 0){
            message(paste("for species", species, sum(occ_data_filtered$ci_width_flag), " rows were filtered as they had confidence intervals that are too large in width"))
        }

        # Provide message for species-years dropped
        if(sum(occ_data_filtered$z_score_flag) > 0){
            message(paste("for species", species, sum(occ_data_filtered$z_score_flag), " rows were filtered as they had z scores that are too large"))
        }

        new_index = occ_data_filtered %>% filter(!ci_width_flag, !z_score_flag) %>% select(-ci_width_flag, -z_score_flag)

        if(nrow(new_index) < min_years_post_filter){
            message(paste("species", species, "has less than", min_years_post_filter , "valid years after filtering. Removing species"))
            species_to_remove = c(species_to_remove, species)
            next
        } else{
            occti_outputs[[species]]$Index <- new_index
        }

    }

    # Remove species that did not have sufficient years of data before or after filtering
    occti_outputs = occti_outputs[!(names(occti_outputs) %in% species_to_remove)]

    return(occti_outputs)

}
