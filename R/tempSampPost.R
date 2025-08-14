#' Extract and Combine Posterior Samples for Multiple Species
#'
#' Iterates over a set of species model outputs, extracts posterior occupancy samples 
#' and associated metadata using \code{\link{combineSamps}}, and returns combined 
#' results for all species in the selection.
#'
#' @param indata Character. Path to the directory containing model output files 
#'   (either \code{.rdata} or \code{.rds} files, depending on how models were saved).
#' @param keep Character vector of species identifiers (file basenames without extensions)
#'   to process. Each entry should match a model file in \code{indata}.
#' @param region Character. Region code for extracting posterior samples (e.g., 
#'   \code{"ENGLAND"}, \code{"WALES"}, or \code{"UK"}). Must match entries in the model 
#'   metadata (\code{out$regions} or \code{out$region_aggs}).
#' @param sample_n Numeric. Target number of posterior draws to retain per species 
#'   (default \code{999}).
#' @param tolerance Numeric. Allowed tolerance between the available number of posterior 
#'   samples and \code{sample_n}. Positive values allow oversampling; negative values 
#'   cause the species to be dropped if fewer than \code{sample_n} draws are available.
#' @param minObs Numeric or \code{NULL}. Minimum number of observations required for 
#'   a species to be included in the output. If \code{NULL}, the check is skipped.
#' @param t0 Numeric. Start year (inclusive) of the analysis window.
#' @param tn Numeric. End year (inclusive) of the analysis window.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{samp_post}}{Data frame of posterior occupancy samples for all species, 
#'     with one row per iteration per species, columns for each year (\code{year_<YYYY>}), 
#'     plus \code{iteration} and \code{species}.}
#'   \item{\code{meta}}{Data frame of combined species-level metadata, including:
#'     \itemize{
#'       \item Observation counts (global and regional)
#'       \item First and last years with observations
#'       \item Largest year gap in observations
#'       \item Rule-of-thumb (ROT) metrics for data quality
#'     }
#'     All metadata columns are suffixed with \code{"_r_<region>"}.
#'   }
#' }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Loops over each species in \code{keep}.
#'   \item Calls \code{\link{combineSamps}} to load the species' model output, 
#'         filter by region and observation count, and extract posterior draws.
#'   \item Binds all posterior samples into one data frame (\code{samp_post}).
#'   \item Binds all metadata rows into one data frame (\code{meta}) and 
#'         appends the region suffix to all column names.
#' }
#'
#' Species that fail to load, lack sufficient observations, or do not meet 
#' \code{sample_n} requirements (within \code{tolerance}) are skipped.
#'
#' @seealso
#' \code{\link{combineSamps}}, \code{\link{applySamp_single}}
#'

tempSampPost <- function(indata,
                         keep,
                         region,
                         sample_n = 999,
                         tolerance = 0,
                         minObs = NULL,
                         t0,
                         tn) {

  run_one <- function(sp) {
    combineSamps(sp,
                 indata = indata,
                 region = region,
                 sample_n = sample_n,
                 tolerance = tolerance,
                 minObs = minObs,
                 scaleObs = scaleObs,
                 t0 = t0,
                 tn = tn)
  }

  outputs <- lapply(keep, run_one)

  samp_post <- do.call("rbind", lapply(outputs, function(x) x[[1]]))
  meta <- do.call("rbind", lapply(outputs, function(x) x[[2]]))

  meta <- data.frame(Species = meta$species,
                     n_obs_global = meta$nRec_glob,
                     n_obs_regional = meta$nRec_reg,
                     min_year_data = meta$first,
                     max_year_data = meta$last,
                     gap_start = 0,
                     gap_end = 0,
                     gap_middle = meta$gap,
                     rot_median = meta$median,
                     rot_P90 = meta$P90,
                     rot_visits_median = meta$visits_median,
                     rot_visits_P90 = meta$visits_P90,
                     rot_prop_list_one = meta$prop_list_one,
                     rot_prop_repeats_grp = meta$prop_repeats_grp,
                     rot_prop_abs = meta$prop_abs,
                     rot_EqualWt = meta$EqualWt,
                     rot_HighSpec = meta$HighSpec)

  colnames(meta) <- paste0(colnames(meta), "_r_", region)

  return(list(samp_post, meta))
}