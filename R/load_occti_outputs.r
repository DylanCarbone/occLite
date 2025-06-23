#' Load OCCTI Occupancy Outputs
#'
#' Loads all OCCTI model output RDS files from a specified directory into a named list.
#'
#' @param results_dir A character string specifying the directory containing OCCTI occupancy output `.rds` files.
#'
#' @return A named list where each element contains the output of one species and is named by the species.
#' 
#' @export
load_occti_outputs <- function(results_dir) {
  # Get full paths to all matching RDS files
  paths <- list.files(results_dir, pattern = "*occupancy_output.rds", full.names = TRUE)
  
  occti_outputs <- list()  # Create an empty list to store the loaded data
  
  for (path in paths) {
    sp_output <- readRDS(path)  # Load each RDS file
    
    species_name <- sp_output$Species  # Use the species name as the list element name
    occti_outputs[[species_name]] <- sp_output  # Store the species output in the list
  }
  
  return(occti_outputs)  # Return the full named list
}