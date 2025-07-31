#' Prepare Species Observation Data for Occupancy Modelling
#'
#' This function reads and preprocesses species observation data by standardising column names, formatting dates,
#' converting grid references to monads, adding spatial coordinates (northing and easting), assigning country region
#' labels (England, Wales, Scotland, or Northern Ireland), and calculating list length per grid square and date.
#'
#' @param data A data frame containing species observation records. Must include columns named `species`, `gridref`, and `date`.
#' @param subset Logical. If \code{TRUE}, the dataset is filtered to retain only species and sites with sufficient data based on \code{min.Recs} and \code{nyr}.
#' @param min.Recs Integer. Minimum number of records per species required for inclusion (only used if \code{subset = TRUE}).
#' @param nyr Integer. Minimum number of years in which a monad must have at least one visit (only used if \code{subset = TRUE}).
#'
#' @details
#' This function performs the following operations:
#' \itemize{
#'   \item Ensures `date` is in \code{Date} format and extracts year and day-of-year.
#'   \item Filters and reformats grid references to monads (1 km precision).
#'   \item Assigns country region labels using an internal reference data frame `monad_region`.
#'   \item Converts monads to British National Grid easting and northing coordinates.
#'   \item Removes records without region information or spatial coordinates.
#'   \item Calls \code{add_dates()} to add derived date-related columns.
#'   \item Calculates list length (`listL`) as the number of unique species per \code{gridref} and \code{date}.
#'   \item Optionally filters data using \code{subset_occ_data()} based on species frequency and visit coverage.
#' }
#'
#' @return A processed data frame containing cleaned and enriched species observations with spatial and temporal fields suitable for occupancy analysis.
#'
#' @import dplyr
#' @importFrom lubridate yday year
#' @importFrom stats na.omit
#' @importFrom BRCmap fmt_gridref reformat_gr OSgrid2GB_EN
#'
#' @export
prep_occ_data = function(data, subset = FALSE, min.Recs = 10, nyr = 2){

# Read and preprocess data
data <- data %>% 
  mutate(date = as.Date(date),
         Year = lubridate::year(date))

# Handle unique grid references to save computation time
unique_gridrefs = unique(data$gridref)

# Remove invalid formats and convert to monad
unique_gridrefs = fmt_gridref(unique_gridrefs)
unique_gridrefs = reformat_gr(unique_gridrefs, prec_out = 1000, pad_gr = FALSE)
unique_gridrefs = unique(unique_gridrefs[!is.na(unique_gridrefs)]) # unique to remove new duplicates

# You have to load the following objects from the BRCmap package
data(datum_vars, package = "BRCmap")
data(helmert_trans_vars, package = "BRCmap")

# Merge with the original data, whilst obtaining the country, and converting to northing and easting
easting_northing_map_df = data.frame(monad = unique_gridrefs) %>%
left_join(monad_region %>% select(monad, region)) %>%
cbind(OSgrid2GB_EN(unique_gridrefs)) %>% 
rename(northing = NORTHING, easting = EASTING)

# Merge with the original data, removing unpaired rows
data = data %>%
left_join(easting_northing_map_df, join_by(gridref == monad)) %>%
filter(!is.na(region))

# Add date variables using occti helper function
data = add_dates(data)

# Calculate list length
data = data %>%
group_by(gridref, date) %>%
mutate(listL = length(unique(species))) %>%
ungroup()

# subset data using subset_occ_data helper function
if (subset){
  data = subset_occ_data(data, min.Recs = min.Recs, nyr = nyr)
}

return(data)
}
