#' @title Pivot occupancy posterior samples to long format
#'
#' @description Converts a wide-format data frame of occupancy model posterior samples 
#' (e.g., from `tempSampPost()`) into a tidy long format, where each row corresponds 
#' to a single prediction for a species-year-iteration combination.
#'
#' @param samp_post A data frame returned by `tempSampPost()` or `applySamp_single()`; 
#' must contain columns named `year_XXXX` for each year, as well as `species` and `iteration`.
#'
#' @return A long-format `data.frame` with columns: `iteration`, `species`, `year`, and `pred`
#'
#' @examples
#' samp_long <- translate_sparta(result$samp_post)
#'
#' @export
translate_sparta <- function(samp_post) {
  samp_post %>%
    tidyr::pivot_longer(
      cols = starts_with("year_"),
      names_to = "year",
      names_prefix = "year_",
      names_transform = list(year = as.integer),
      values_to = "pred"
    )
}