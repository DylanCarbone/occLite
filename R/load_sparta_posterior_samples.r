#' Load and process SPARTA posterior samples
#'
#' Loads a saved \code{.RData} file containing SPARTA posterior samples, 
#' extracts the first element of the loaded object (assumed to be a data frame), 
#' and converts the \code{species} column to lower case with spaces replaced 
#' by underscores.
#'
#' @param posterior_samples_file Character string. Path to the \code{.RData} file 
#'   containing the posterior samples object (expected to be named \code{comb} in 
#'   the workspace).
#'
#' @return A \code{data.frame} or \code{tibble} containing the processed posterior 
#'   samples for the first element of \code{comb}.
#'
#' @details
#' This function assumes the \code{.RData} file contains an object named \code{comb}, 
#' which is a list whose first element contains a column \code{species}. The 
#' \code{species} values are transformed to lower case and spaces are replaced with 
#' underscores for consistency.
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' posterior_df <- load_sparta_posterior_samples("data/sparta_samples.RData")
#' }
#'
#' @export
#' 
load_sparta_posterior_samples <- function(posterior_samples_file) {
  
  # Load the RData file and retrieve the object it contains
  load(posterior_samples_file)
  
  # Process the first element of the list
  comb[[1]] %>%
    dplyr::mutate(
      species = tolower(gsub(" ", "_", species))
    )
}
