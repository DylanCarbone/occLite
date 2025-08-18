#' Load OCCTI Occupancy Outputs
#'
#' Loads all OCCTI model output RDS files from a specified directory into a named list.
#'
#' @param results_dir A character string specifying the directory containing OCCTI occupancy output `.rds` files.
#'
#' @return A named list where each element contains the output of one species and is named by the species.
#' 
#' @export
load_occti_outputs <- function(results_dir, pattern = "*occupancy_output.rds") {

  # Get full paths to all matching RDS files
  paths <- list.files(results_dir, pattern = pattern, full.names = TRUE)

  if(length(paths) == 0){
    stop("No files found. Are you sure you specified the correct directory and pattern?")
  }

  occti_outputs <- list()  # Create an empty list to store the loaded data
  
  for (path in paths) {
    sp_output <- readRDS(path)  # Load each RDS file

    species_name <- as.character(sp_output$species)

    if(!is.null(occti_outputs[[species_name]])){
      stop("There are duplicate species files. Did you subset your pattern by region?")
    }
      
    occti_outputs[[species_name]] <- sp_output  # Store the species output in the list
  }
  
  return(occti_outputs)  # Return the full named list
}