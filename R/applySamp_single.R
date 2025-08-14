#' Apply Posterior Sampling for a Single Region
#'
#' Runs the posterior sampling process for all species models in a given region, 
#' applying observation count filters and optionally restricting to a user-specified 
#' set of species. This function is a wrapper around \code{tempSampPost()} 
#' for a single region, reading model outputs from \code{.rdata} files and returning 
#' sampled posterior draws and associated metadata.
#'
#' @param modPath Character. Path to the directory containing model output files 
#'   in \code{.rdata} format. Each file should correspond to a single species.
#' @param region Character or \code{NULL}. Region code for filtering results. 
#'   If \code{NULL}, all models are processed without regional subsetting. 
#'   Should match the region names used in the model metadata.
#' @param nSamps Numeric. Number of posterior samples to draw per species 
#'   (default \code{999}).
#' @param minObs Numeric. Minimum number of observations required per species 
#'   within the time window \code{t0}–\code{tn} for inclusion in sampling 
#'   (default \code{50}).
#' @param tolerance Numeric. Allowed tolerance when matching the available 
#'   number of iterations to \code{nSamps}. Positive values allow slight over-sampling, 
#'   negative values will trigger an error if fewer than \code{nSamps} iterations are available 
#'   (default \code{0}).
#' @param write Logical. If \code{TRUE}, saves the combined sampling output 
#'   (posterior samples and metadata) to \code{outPath}.
#' @param outPath Character. File path to save results when \code{write = TRUE}.
#'   The object is saved as an R binary file (\code{.RData}) containing a list 
#'   with two elements: sampled posteriors and metadata.
#' @param speciesToKeep Character vector (comma-separated string with no spaces) or \code{NA}.
#'   Optional list of species names to include. Matching is case-insensitive. 
#'   Species not found in \code{modPath} will trigger a warning.
#' @param t0 Numeric. First year of the analysis window (inclusive).
#' @param tn Numeric. Last year of the analysis window (inclusive).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{samp_post}}{Data frame of posterior samples for each species, 
#'   with columns for species name, iteration number, and posterior estimates per year.}
#'   \item{\code{meta}}{Data frame of metadata for each species, including 
#'   record counts, first and last years with data, and rule-of-thumb metrics.}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Identifies available \code{.rdata} model files in \code{modPath}.
#'   \item Optionally filters to a given set of species.
#'   \item Runs \code{tempSampPost()} to extract and sample posterior draws.
#'   \item Optionally writes the results to disk.
#' }
#' The filtering step uses \code{minObs} to drop species with insufficient 
#' observations in the specified time window (\code{t0}–\code{tn}). 
#' Sampling is performed if more posterior draws are available than \code{nSamps}, 
#' subject to the \code{tolerance} setting.
#'
#' @seealso \code{\link{tempSampPost}}, \code{\link{combineSamps}}
#'
#' @export

applySamp_single <- function(modPath, region = NULL,
                             nSamps = 999, minObs = 50,
                             tolerance = 0,
                             write = FALSE,
                             outPath = NULL, speciesToKeep = NA,
                             t0, tn) {

  if(is.null(region)){
    message("You have specified region as NULL. To change this, please modify the 'region' parameter")
  }

  # ----- Identify available model files -----
  modFiles_rdata <- list.files(modPath, pattern = ".rdata")
  
  if (length(modFiles_rdata) == 0){
    stop("Model files must be .rdata")
  }

  keep = gsub(".rdata", "", modFiles_rdata)
  
  # ----- Filter to speciesToKeep if provided -----
  if (!is.na(speciesToKeep)) {
    speciesToKeep_vec <- unlist(strsplit(speciesToKeep, ","))
    notFound <- speciesToKeep_vec[!tolower(speciesToKeep_vec) %in% tolower(keep)]
    
    if (length(notFound) > 0) {
      warning("Some species in 'speciesToKeep' not found: ", paste(notFound, collapse = ", "))
    }
    
    keep <- keep[tolower(keep) %in% tolower(speciesToKeep_vec)]
  }
  
  # ----- Run sampling or extract means -----
  out <- tempSampPost(indata = modPath,
                      keep = keep,
                      region = region,
                      sample_n = nSamps,
                      tolerance = tolerance,
                      minObs = minObs,
                      t0 = t0,
                      tn = tn)

  samp_post <- out[[1]]
  meta <- out[[2]]
  
  samp_post$species <- tolower(samp_post$species)
  meta[, 1] <- tolower(meta[, 1])
  
  if (write && !is.null(outPath)) {
    comb <- list(samp_post, meta)
    save(comb, file = outPath)
  }
  
  return(out)
}
