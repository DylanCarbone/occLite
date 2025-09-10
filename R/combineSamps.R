#' Load and Sample Posterior Distributions for a Single Species
#'
#' Loads model output for a given species from an \code{.rdata} file, extracts the 
#' posterior occupancy samples for a specified region (or globally), filters the data 
#' by a specified time window and minimum observation threshold, and optionally 
#' downsamples the posterior draws to a target number.
#'
#' @param species Character. Species name or identifier corresponding to the 
#'   model output file \code{<species>.rdata} in \code{indata}.
#' @param indata Character. Path to the directory containing \code{.rdata} model output files.
#' @param region Character or \code{NULL}. Region code for extracting posterior 
#'   samples (e.g., \code{"ENGLAND"}, \code{"WALES"}, etc.). If \code{NULL}, 
#'   global samples are extracted.
#' @param sample_n Numeric. Target number of posterior draws to retain.
#' @param tolerance Numeric. Allowed tolerance between the number of available 
#'   posterior samples and \code{sample_n}. Positive values allow slight over-sampling; 
#'   negative values will cause the function to stop if fewer than \code{sample_n} 
#'   samples are available.
#' @param minObs Numeric or \code{NULL}. Minimum number of observations required 
#'   within the analysis window (\code{t0}–\code{tn}) for the species to be included. 
#'   If \code{NULL}, the check is skipped.
#' @param t0 Numeric. Start year of the analysis window (inclusive).
#' @param tn Numeric. End year of the analysis window (inclusive).
#'
#' @return A list of length 2:
#' \describe{
#'   \item{\code{raw_occ}}{Data frame containing posterior occupancy samples for 
#'     the species, with one column per year, plus \code{iteration} and \code{species} 
#'     columns.}
#'   \item{\code{meta}}{Data frame of summary metadata for the species, including 
#'     record counts (\code{nRec_glob}, \code{nRec_reg}), first and last observation 
#'     years, largest observation gap, analysis window, and placeholder columns 
#'     for rule-of-thumb metrics.}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Loads the model output (\code{.rdata}) for the given species.
#'   \item Extracts observation data and filters to the time window (\code{t0}–\code{tn}).
#'   \item If \code{region} is specified:
#'     \itemize{
#'       \item Uses \code{out$regions} if it is a direct region.
#'       \item Uses \code{out$region_aggs} if it is an aggregated region (combining multiple subregions).
#'     }
#'   \item Computes the number of observations (\code{nRec}) and drops the species if 
#'     \code{nRec < minObs}.
#'   \item Checks that the requested posterior variable exists in \code{out$BUGSoutput$sims.list}.
#'   \item Samples or truncates posterior draws to \code{sample_n}, within the given \code{tolerance}.
#'   \item Constructs a metadata row summarising observation counts, years, and gaps.
#' }
#'
#' The posterior samples are returned in wide format, with one column per model year 
#' prefixed with \code{"year_"}.
#'
#' @seealso
#' \code{\link{applySamp_single}}, \code{\link{tempSampPost}}
#'
#' @examples
#' \dontrun{
#' res <- combineSamps(
#'   species = "apis_mellifera",
#'   indata = "models/",
#'   region = "ENGLAND",
#'   sample_n = 500,
#'   tolerance = 10,
#'   minObs = 30,
#'   t0 = 2000,
#'   tn = 2020
#' )
#' }
#'

combineSamps <- function(species, 
                           indata, 
                           region, 
                           sample_n, 
                           tolerance,
                           minObs,
                           t0, 
                           tn) {

  region_psi_fs <- paste0("psi.fs", ifelse(!is.null(region), paste0(".r_", region), ""))

  filepath <- file.path(indata, paste0(species, ".rdata"))
  out <- NULL

  load_attempt <- try(load(filepath), silent = TRUE)

  if (inherits(load_attempt, "try-error") || is.null(out$model)) {
    message("Model failed to load or has no valid model object: ", species)
    return(NULL)
  }

  dat <- out$model$data()
  dat_df <- data.frame(year = dat$Year, rec = dat$y, site = dat$Site)

  if (is.null(region)) {
    dat_glob <- dat_df[dat_df$year >= (t0 - (out$min_year - 1)) & dat_df$year <= (tn - (out$min_year - 1)), ]
    nRec_glob <- sum(dat_glob$rec)
    nRec <- nRec_glob
    datm <- dat_glob
  } else {
    if (region %in% out$regions) {
      region_site <- dat[[paste0("r_", region)]][dat$Site]
    } else {
      region_aggs <- unlist(out$region_aggs[region])
      region_site <- rowSums(sapply(region_aggs, function(x) dat[[paste0("r_", x)]][dat$Site]))
    }

    dat_reg <- cbind(dat_df, region_site)
    dat_reg <- dat_reg[dat_reg$region_site == 1 & dat_reg$year >= (t0 - (out$min_year - 1)) & dat_reg$year <= (tn - (out$min_year - 1)), ]
    nRec_reg <- sum(dat_reg$rec)
    nRec <- nRec_reg
    datm <- dat_reg
  }
  if (!is.null(minObs) && nRec < minObs) {
    message("Dropped (not enough observations): ", species)
    return(NULL)
  }

  if (!region_psi_fs %in% names(out$BUGSoutput$sims.list)) {
    message("Posterior region variable not found: ", region_psi_fs)
    return(NULL)
  }

  raw_occ <- as.data.frame(out$BUGSoutput$sims.list[[region_psi_fs]])
  if (nrow(raw_occ) == 0) {
    message("Empty posterior sample matrix for: ", species)
    return(NULL)
  }

  diff <- nrow(raw_occ) - sample_n

  n <- nrow(raw_occ)
  if (n >= sample_n) {
    raw_occ <- raw_occ[sample(seq_len(n), sample_n), , drop = FALSE]
  } else {
    short_by <- sample_n - n
    if (short_by > tolerance) {
      message("Not enough posterior samples for ", species,
              " (have ", n, ", need ", sample_n, ", tolerance ", tolerance, ").")
      return(NULL)
    }
  }

  first <- min(datm$year[datm$rec == 1]) + (out$min_year - 1)
  last <- max(datm$year[datm$rec == 1]) + (out$min_year - 1)

  yrs <- sort(unique(datm$year[datm$rec == 1]))
  gap <- if (length(yrs) > 1) max(diff(yrs)) else 1   # largest run of missing years
  
  # --- Extract rule-of-thumb metrics if present ---
  rot <- NULL
  if (!is.null(attr(out, "metadata")$analysis$spp_Metrics)) {
    rot <- as.data.frame(attr(out, "metadata")$analysis$spp_Metrics)

    # EqualWt and HighSpec decision trees (see https://www.biorxiv.org/content/10.1101/813626v1.full)
    rot$EqualWt <- ifelse(rot$prop_abs >= 0.990, rot$P90 >= 3.1, rot$P90 >= 6.7)
    rot$HighSpec <- ifelse(rot$prop_abs >= 0.958, rot$P90 >= 9.5, rot$P90 >= 29)
  }else{
    rot <- data.frame(median = NA, P90 = NA, visits_median = NA, visits_P90 = NA, prop_list_one = NA, prop_repeats_grp = NA, prop_abs = NA, EqualWt = NA, HighSpec = NA)
  }

  # construct metadata table
  meta <- data.frame(species = species,
                     nRec_glob = ifelse(exists("nRec_glob"), nRec_glob, NA),
                     nRec_reg = ifelse(exists("nRec_reg"), nRec_reg, NA),
                     first = first, last = last,
                     gap = gap,
                     firstMod = t0, lastMod = tn,
                     rot,
                     row.names = NULL)

  # --- Rename posterior columns consistently ---
  colnames(raw_occ) <- paste0("year_", out$min_year:out$max_year)
  
  # --- Always add iteration and species as named columns ---
  raw_occ$iteration <- seq_len(nrow(raw_occ))
  raw_occ$species   <- species
  
  # Restrict posterior samples to t0–tn
  keep_years <- paste0("year_", t0:tn)
  keep_years <- keep_years[keep_years %in% names(raw_occ)]
  raw_occ <- raw_occ[, c(keep_years, "iteration", "species"), drop = FALSE]

  return(list(raw_occ = raw_occ, meta = meta))
}
