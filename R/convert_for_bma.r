#' Convert LOESS Summary List to BMA-Compatible Data Frame
#'
#' This function converts the LOESS summary outputs from `smooth_occti_outputs()` into a format
#' suitable for Bayesian Model Averaging (BMA) input. It binds all species data into a single
#' data frame and renames key columns to match expected BMA input structure.
#'
#' @param loess_summaries A list object returned by `smooth_occti_outputs()`, 
#'   which must contain a `loess_summaries` element (a named list of data frames per species).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{species}{Species name from the list names}
#'     \item{year}{Year of the estimate}
#'     \item{index}{Smoothed occupancy index (mean)}
#'     \item{se}{Standard error of the smoothed index}
#'   }
#'
#' @import gtools
#' @importFrom dplyr bind_rows rename select
#' @export
convert_for_bma <- function(loess_summaries) {
  
  # Combine all LOESS summaries into one data frame, keeping species names as an identifier
  loess_summaries_df <- bind_rows(loess_summaries$loess_summaries, .id = "species")

  # Rename columns to match BMA input expectations and select relevant ones
  loess_summaries_df_named <- loess_summaries_df %>%
    rename(
      year = Year,
      index = psiA_loess_mean,
      se = psiA_loess_se
    ) %>%
    select(species, year, index, se)

  return(loess_summaries_df_named)
}
