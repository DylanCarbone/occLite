#' applySamp_single - Applies the sampling process for a single taxonomic group and filter set.
#'
#' @param modPath Character. Base path to the model outputs.
#' @param metaPath Character. Path to the metadata.
#' @param group Character. Taxonomic group (e.g. "Bees").
#' @param region Character. One of "UK", "GB", "ENGLAND", "WALES", "SCOTLAND", or "NORTHERN_IRELAND".
#' @param indicator Character. One of "priority", "pollinators", or "all".
#' @param nSamps Numeric. Number of posterior samples to draw.
#' @param minObs Numeric. Minimum number of observations per species.
#' @param scaleObs Character. Either "global" or "region".
#' @param sample Logical. If TRUE, samples the posterior; if FALSE, just extracts mean estimates.
#' @param write Logical. If TRUE, saves results to `outPath`.
#' @param outPath Character. Path to save output if `write = TRUE`.
#' @param speciesToKeep Character vector of species to keep (comma-separated string or NA).
#' @param drop Logical. Whether to drop species based on exclusion advice.
#' @param t0, tn Numeric. Start and end years for analysis.
#' @param clipBy Character. One of "species" or "group".
#' @param parallel Logical. Whether to use parallel processing.
#'
#' @return A list with elements: samp_post, meta, indicator, group_name, region, clipBy, minObs.
#'
#' @export

applySamp_single <- function(modPath, metaPath, group, region,
                             indicator = "all", nSamps = 999, minObs = 50,
                             scaleObs = "global", sample = TRUE, write = FALSE,
                             outPath = NULL, speciesToKeep = NA,
                             drop = TRUE, t0, tn, clipBy = "group", parallel = TRUE) {
  
  data(speciesInfo)
  
  keep <- NULL
  keep_iter <- NULL
  
  # ----- Priority species filtering -----
  if (indicator == "priority") {
    if (!region %in% c("WAL", "WALES", "SCO", "SCOTLAND", "ENG", "ENGLAND", "NIR", "NORTHERN_IRELAND", "GB", "UK")) {
      stop("Priority species region not recognised")
    }
    
    pr_region <- switch(region,
                        WAL = "WALES", WALES = "WALES",
                        ENG = "ENGLAND", ENGLAND = "ENGLAND",
                        SCO = "SCOTLAND", SCOTLAND = "SCOTLAND",
                        NIR = "NORTHERN.IRELAND", NORTHERN_IRELAND = "NORTHERN.IRELAND",
                        GB = c("WALES", "ENGLAND", "SCOTLAND"),
                        UK = c("WALES", "ENGLAND", "SCOTLAND", "NORTHERN.IRELAND"))
    
    speciesInfo_group <- speciesInfo[speciesInfo$Group %in% group, ]
    
    keepInds <- sapply(1:nrow(speciesInfo_group), function(i) any(speciesInfo_group[i, pr_region] == "Y"))
    
    keep <- c(as.character(speciesInfo_group$Species[keepInds]),
              as.character(speciesInfo_group$concept[keepInds]))
    
    keep <- unique(keep)
  }
  
  # ----- Identify available model files -----
  modFiles_rdata <- list.files(modPath, pattern = ".rdata")
  modFiles_rds   <- list.files(modPath, pattern = ".rds")
  
  if (length(modFiles_rdata) == 0 && length(modFiles_rds) == 0) {
    stop("Model files must be either .rds or .rdata")
  }
  
  filetype <- if (length(modFiles_rdata) > 0) "rdata" else "rds"
  modFiles <- if (filetype == "rdata") modFiles_rdata else modFiles_rds
  
  # keep_iter <- gsub(paste0("\\.", filetype), "", modFiles)
  
  # ----- Derive species names for non-priority cases -----
  if (indicator != "priority") {
    # first_spp <- keep_iter[[1]]
    
    # if (grepl("_[0-9]+$", first_spp)) {
    #   keep <- gsub("(.*)_\\w+", "\\1", keep_iter)
    #   keep <- unique(keep)
    # } else {
    #   keep <- keep_iter
    #   keep_iter <- NULL
    # }

    # The above is not necessary because we are not chaining
    keep = gsub(paste0(".", filetype), "", modFiles)
  }
  
  # ----- Filter to speciesToKeep if provided -----
  if (!is.na(speciesToKeep)) {
    speciesToKeep_vec <- unlist(strsplit(speciesToKeep, ","))
    notFound <- speciesToKeep_vec[!tolower(speciesToKeep_vec) %in% tolower(keep)]
    
    if (length(notFound) > 0) {
      warning("Some species in 'speciesToKeep' not found: ", paste(notFound, collapse = ", "))
    }
    
    keep <- keep[tolower(keep) %in% tolower(speciesToKeep_vec)]
  }
  
  # ----- Drop excluded species -----
  if (drop) {
    drop_ids <- which(!is.na(speciesInfo$Reason_not_included) & speciesInfo$Reason_not_included != "Didn't meet criteria")
    drop_list <- c(as.character(speciesInfo$Species[drop_ids]), as.character(speciesInfo$concept[drop_ids]))
    keep <- keep[!keep %in% drop_list]
  }
  
  # ----- Run sampling or extract means -----
  if (sample) {
    out <- tempSampPost(indata = modPath,
                        keep = keep,
                        keep_iter = NULL,
                        output_path = NULL,
                        region = region,
                        sample_n = nSamps,
                        group_name = group,
                        combined_output = TRUE,
                        write = FALSE,
                        minObs = minObs,
                        scaleObs = scaleObs,
                        t0 = t0,
                        tn = tn,
                        parallel = parallel,
                        filetype = filetype)
  } else {
    out <- getA(indata = modPath,
                keep = keep,
                REGION_IN_Q = paste0("a_", region),
                group_name = group,
                combined_output = TRUE,
                write = FALSE,
                minObs = minObs,
                t0 = t0,
                tn = tn,
                parallel = parallel)
  }
  
  samp_post <- out[[1]]
  meta <- out[[2]]
  
  samp_post$species <- tolower(samp_post$species)
  meta[, 1] <- tolower(meta[, 1])
  
  if (write && !is.null(outPath)) {
    comb <- list(samp_post, meta)
    save(comb, file = file.path(outPath, paste0(group, "_", indicator, "_", region, "_samp.rdata")))
  }
  
  return(list(
    samp_post = samp_post,
    meta = meta,
    indicator = indicator,
    group_name = group,
    region = region,
    clipBy = clipBy,
    minObs = minObs
  ))
}
