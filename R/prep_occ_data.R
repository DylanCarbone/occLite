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

if(!all(c("date", "species", "gridref") %in% colnames(data))){
  stop("Your dataframe must include columns named date, species, and gridref")
}

# Start with unique gridref
new_gridref <- distinct(data, gridref)

# --- Step 1: formatting ---
new_gridref <- new_gridref %>%
  mutate(grid_refs_post_formatting = fmt_gridref(gridref))

n_fail1 <- sum(is.na(new_gridref$grid_refs_post_formatting))

message(n_fail1, " unique grid references have invalid formats")

# --- Step 2: resolution standardisation ---
new_gridref <- new_gridref %>%
  filter(!is.na(grid_refs_post_formatting)) %>%
  mutate(grid_refs_post_resolution_standardisation =
           reformat_gr(grid_refs_post_formatting,
                       prec_out = 1000, pad_gr = FALSE))

n_fail2 <- sum(is.na(new_gridref$grid_refs_post_resolution_standardisation))

message(n_fail2, " unique grid references have a resolution greater than monad resolution")

# --- Step 3: region join ---
new_gridref <- new_gridref %>%
  filter(!is.na(grid_refs_post_resolution_standardisation)) %>%
  left_join(monad_region %>% select(monad, region),
            join_by(grid_refs_post_resolution_standardisation == monad))

n_fail3 <- sum(is.na(new_gridref$region))

message(n_fail3, " unique grid references could not be mapped to a region, as they fall outside of the UK")

# --- Step 4: UK coordinate conversion ---

# Filter grid references that have not been mapped to a region
new_gridref <- new_gridref %>%
    filter(!is.na(region))

# You have to load the following objects from the BRCmap package
data(datum_vars, package = "BRCmap")
data(helmert_trans_vars, package = "BRCmap")

new_gridref <- cbind(new_gridref, OSgrid2GB_EN(new_gridref$grid_refs_post_resolution_standardisation)) %>% 
    rename(northing = NORTHING, easting = EASTING) %>%
    select(-grid_refs_post_formatting)

# --- Final merge ---
data <- data %>%
  left_join(new_gridref, by = "gridref") %>%
  filter(!is.na(region))  %>% # This filters out all invalid grid references from the previous steps
    select(-gridref)  %>%
    rename(gridref = grid_refs_post_resolution_standardisation)

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
