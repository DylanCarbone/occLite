#' Filter posterior samples by RoT rules and report removals
#'
#' @param posterior_samples A list of length 2 containing posterior samples and sparta run metadata
#'
#' @return A list with elements:
#'   - posterior_samples: filtered samples
#'   - metadata: filtered metadata (only kept species)
#'
#' @import dplyr
#' @export 
#' 

filter_rot <- function(posterior_samples) {

  samples  <- posterior_samples$posterior_samples
  metadata <- posterior_samples$metadata

  # Stage 1: Filter by rot_equalwt_r_uk
  n_start <- nrow(metadata)
  metadata <- metadata %>% filter(rot_equalwt_r_uk == TRUE)
  message(sprintf("%d species removed as rot_equalwt_r_uk was NA or FALSE", n_start - nrow(metadata)))

  if(nrow(metadata) == 0){
    stop("Your metadata has zero rows after filtering for TRUE rot_equalwt_r_uk values. Filtering by rule of thumb is not possible")
  }

  # Stage 2: Apply pass rule
  n_start <- nrow(metadata)
  metadata <- metadata %>%
    mutate(pass = ifelse(rot_prop_abs_r_uk >= 0.990, 
                         rot_p90_r_uk >= 3.1, 
                         rot_p90_r_uk >= 6.7)) %>%
    filter(pass)
  
  message(sprintf("%d species removed as they did not pass rule of thumb threshold", n_start - nrow(metadata)))

  # Filter samples
  samples  <- samples %>% filter(species_r_uk %in% metadata$species_r_uk)

  return(list(posterior_samples = samples, metadata = metadata))
}