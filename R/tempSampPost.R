#' tempSampPost - updated version to handle flat filenames (e.g. tik_123.rdata)
#'
#' @param indata Path to the directory containing .rdata model files
#' @param keep Character vector of file basenames without extension (e.g. "tik_123")
#' @param keep_iter NULL unless using chained models
#' @param output_path Where to save combined posterior samples (if write = TRUE)
#' @param region Region label (used in metadata column names)
#' @param sample_n Number of posterior samples to extract (default = 999)
#' @param tolerance Allowed deviation from sample_n
#' @param group_name Group name (used in output filename)
#' @param combined_output Logical, whether to merge posterior samples across species
#' @param max_year_model Optionally filter model years
#' @param min_year_model Optionally filter model years
#' @param write Logical, whether to write output to disk
#' @param minObs Filter species with fewer observations (optional)
#' @param scaleObs Either "global" or "regional"
#' @param t0 Start year for summary
#' @param tn End year for summary
#' @param parallel Logical, whether to use parallel execution
#' @param n.cores Optional override for number of cores
#' @param filetype Either "rdata" or "rds"
#'
#' @return List with posterior samples and metadata
#' @export

tempSampPost <- function(indata,
                         keep,
                         keep_iter = NULL,
                         output_path,
                         region,
                         sample_n = 999,
                         tolerance = 0,
                         group_name = "",
                         combined_output = TRUE,
                         max_year_model = NULL,
                         min_year_model = NULL,
                         write = FALSE,
                         minObs = NULL,
                         scaleObs = "global",
                         t0,
                         tn,
                         parallel = FALSE,
                         n.cores = NULL,
                         filetype = "rdata") {

  if (parallel && is.null(n.cores)) {
    n.cores <- parallel::detectCores() - 1
  }

  iter <- NULL
  if (!is.null(keep_iter)) {
    findIteration <- function(list_of_file_names) {
      list_of_file_names <- gsub('_[[:digit:]]{1}$', '', list_of_file_names)
      iterations <- regmatches(list_of_file_names, regexpr('[[:digit:]]+$', list_of_file_names))
      return(c(min(as.numeric(iterations)), max(as.numeric(iterations))))
    }
    iter <- findIteration(keep_iter)
  }

  run_one <- function(sp) {
    combineSamps(sp,
                 indata = indata,
                 keep_iter = keep_iter,
                 region = region,
                 sample_n = sample_n,
                 tolerance = tolerance,
                 combined_output = combined_output,
                 minObs = minObs,
                 scaleObs = scaleObs,
                 t0 = t0,
                 tn = tn,
                 filetype = filetype,
                 iter = iter)
  }

  outputs <- if (parallel) {
    parallel::mclapply(keep, run_one, mc.cores = n.cores)
  } else {
    lapply(keep, run_one)
  }

  samp_post <- do.call("rbind", lapply(outputs, function(x) x[[1]]))
  meta <- do.call("rbind", lapply(outputs, function(x) x[[2]]))

  meta <- data.frame(Species = meta$species,
                     n_obs_global = meta$nRec_glob,
                     n_obs_regional = meta$nRec_reg,
                     min_year_data = meta$first,
                     max_year_data = meta$last,
                     min_year_model = meta$firstMod,
                     max_year_model = meta$lastMod,
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

  if (write) {
    if (is.null(output_path)) stop("Must provide output_path when write = TRUE")
    save(samp_post, file = file.path(output_path, paste0(group_name, "_all_spp_sample_", sample_n, "_post_", region, ".rdata")))
  }

  return(list(samp_post, meta))
}