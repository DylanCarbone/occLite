#' Subset Occupancy Data by Species and Site Sampling Frequency
#'
#' Filters an occupancy dataset to include only species with a minimum number of records,
#' and sites that have been sampled in at least a specified number of distinct years.
#'
#' This is useful for cleaning raw species observation data before fitting occupancy models.
#'
#' @param data A `data.frame` or tibble containing the columns `species`, `gridref`, and `Year`.
#' @param min.Recs An integer specifying the minimum number of records required for a species to be retained. Default is 10.
#' @param nyr An integer specifying the minimum number of distinct years a site must be sampled to be retained. Default is 2.
#'
#' @return A filtered `data.frame` containing only rows meeting both species and site sampling criteria.
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' cleaned_data <- subset_occ_data(data, min.Recs = 50, nyr = 3)
#' }
#'
#' @export
subset_occ_data <- function(data, min.Recs = 10, nyr = 2) {
  # Store original numbers for reporting
  original_n <- nrow(data)
  original_species <- length(unique(data$species))
  original_sites <- length(unique(data$gridref))
  
  # Summarise species
  speciesSummary <- data %>%
    group_by(species) %>%
    summarise(
      nuSiteID = n_distinct(gridref),
      nuRecs = n(),
      .groups = "drop"
    )
  
  # Filter by species
  allSpecies <- sort(speciesSummary$species[speciesSummary$nuRecs > min.Recs])
  data_species_filtered <- data %>% filter(species %in% allSpecies)
  species_filtered_n <- nrow(data_species_filtered)
  species_removed <- original_species - length(unique(data_species_filtered$species))
  
  cat("Filtered by species with >", min.Recs, "records:\n")
  cat(species_removed, "species removed\n")
  cat(original_n - species_filtered_n, "rows removed\n\n")
  
  # Filter by site sampling frequency
  sites_to_include <- data %>%
    distinct(gridref, Year) %>%
    group_by(gridref) %>%
    summarise(n_years_sampled = n(), .groups = "drop") %>%
    filter(n_years_sampled >= nyr) %>%
    pull(gridref)
  
  final_data <- data_species_filtered %>%
    filter(gridref %in% sites_to_include)
  
  final_n <- nrow(final_data)
  final_species <- length(unique(final_data$species))
  final_sites <- length(unique(final_data$gridref))
  sites_removed <- original_sites - final_sites
  
  cat("Filtered by sites sampled over >=", nyr, "years:\n")
  cat(sites_removed, "sites removed\n")
  cat(species_filtered_n - final_n, "rows removed\n\n")
  
  return(final_data)
}
