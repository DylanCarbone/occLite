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

filter_rot <- function(posterior_samples, region = "uk") {

  samples  <- posterior_samples$posterior_samples
  metadata <- posterior_samples$metadata

  # Stage 1: Filter by rot_equalwt_r_<region>
  n_start <- nrow(metadata)

  # Make sure region is lowercase
  region = tolower(region)

  metadata <- metadata %>% filter(!!sym(paste0("rot_equalwt_r_", region)) == TRUE)
  message(sprintf("%d species removed as rot_equalwt_r_%s was NA or FALSE", n_start - nrow(metadata), region))

  if(nrow(metadata) == 0){
    stop(paste("Your metadata has zero rows after filtering for TRUE", paste0("rot_equalwt_r_", region), "values. Filtering by rule of thumb is not possible"))
  }
  
  # Stage 2: Apply pass rule
  n_start <- nrow(metadata)
  metadata <- metadata %>%
    mutate(pass = ifelse(!!sym(paste0("rot_prop_abs_r_", region)) >= 0.990, 
                         !!sym(paste0("rot_p90_r_", region)) >= 3.1, 
                         !!sym(paste0("rot_p90_r_", region)) >= 6.7)) %>%
    filter(pass)
  
  message(sprintf("%d species removed as they did not pass rule of thumb threshold", n_start - nrow(metadata)))

  # Filter samples
  samples  <- samples %>% filter(species %in% metadata[[paste0("species_r_", region)]])

  return(list(posterior_samples = samples, metadata = metadata))
}