#' combineSamps
#'
#' This function loads and samples posterior distributions from model outputs for a given species.
#' It is based on the original `combineSamps` function used in `tempSampPost` from the wrappeR package:
#' https://github.com/BiologicalRecordsCentre/wrappeR
#'
#' @export
combineSamps <- function(species, 
                           indata, 
                           keep_iter, 
                           region, 
                           sample_n, 
                           tolerance,
                           combined_output,
                           minObs, 
                           scaleObs,
                           t0, 
                           tn, 
                           filetype,
                           iter) {

  REGION_IN_Q <- paste0("psi.fs.r_", region)

  load_rdata <- function(fileName) {
    load(fileName)
    get(ls()[ls() != "fileName"])
  }

  filepath <- file.path(indata, paste0(species, ".", filetype))
  out_dat <- NULL
  out_meta <- NULL

  if (!is.null(keep_iter)) {
    min_iter <- iter[1]
    max_iter <- iter[2]
    out_dat <- try(load_rdata(file.path(indata, paste0(species, "_", max_iter, "_1.rdata"))), silent = TRUE)
    out_meta <- try(load_rdata(file.path(indata, paste0(species, "_", min_iter, "_1.rdata"))), silent = TRUE)
    if (!inherits(out_dat, "try-error") && is.null(out_dat$model)) out_dat <- out_dat$out
  } else {
    if (filetype == "rds") {
      out_dat <- try(readRDS(filepath), silent = TRUE)
      out_meta <- out_dat
    } else if (filetype == "rdata") {
      out_dat <- try(load_rdata(filepath), silent = TRUE)
      out_meta <- out_dat
    }
  }

  if (inherits(out_dat, "try-error") || is.null(out_dat$model)) {
    message("Model failed to load or has no valid model object: ", species)
    return(NULL)
  }

  dat <- out_meta$model$data()
  dat_df <- data.frame(year = dat$Year, rec = dat$y, site = dat$Site)

  if (scaleObs == "global") {
    dat_glob <- dat_df[dat_df$year >= (t0 - (out_meta$min_year - 1)) & dat_df$year <= (tn - (out_meta$min_year - 1)), ]
    nRec_glob <- sum(dat_glob$rec)
    nRec <- nRec_glob
    datm <- dat_glob
  } else {
    if (region %in% out_meta$regions) {
      region_site <- dat[[paste0("r_", region)]][dat$Site]
    } else {
      region_aggs <- unlist(out_meta$region_aggs[region])
      region_site <- rowSums(sapply(region_aggs, function(x) dat[[paste0("r_", x)]][dat$Site]))
    }
    dat_reg <- cbind(dat_df, region_site)
    dat_reg <- dat_reg[dat_reg$region_site == 1 & dat_reg$year >= (t0 - (out_meta$min_year - 1)) & dat_reg$year <= (tn - (out_meta$min_year - 1)), ]
    nRec_reg <- sum(dat_reg$rec)
    nRec <- nRec_reg
    datm <- dat_reg
  }

  if (!is.null(minObs) && nRec < minObs) {
    message("Dropped (not enough observations): ", species)
    return(NULL)
  }

  if (!REGION_IN_Q %in% names(out_dat$BUGSoutput$sims.list)) {
    message("Posterior region variable not found: ", REGION_IN_Q)
    return(NULL)
  }

  raw_occ <- as.data.frame(out_dat$BUGSoutput$sims.list[[REGION_IN_Q]])
  if (nrow(raw_occ) == 0) {
    message("Empty posterior sample matrix for: ", species)
    return(NULL)
  }

  diff <- nrow(raw_occ) - sample_n

  if (diff < -tolerance) {
    message("Not enough posterior samples for ", species)
    return(NULL)
  } else if (diff > tolerance) {
    raw_occ <- raw_occ[sample(seq_len(nrow(raw_occ)), sample_n), ]
  } else {
    raw_occ <- raw_occ[seq_len(min(nrow(raw_occ), sample_n)), ]
  }

  colnames(raw_occ) <- paste0("year_", out_meta$min_year:out_meta$max_year)
  raw_occ$iteration <- seq_len(nrow(raw_occ))
  raw_occ$species <- species

  first <- min(datm$year[datm$rec == 1]) + (out_meta$min_year - 1)
  last <- max(datm$year[datm$rec == 1]) + (out_meta$min_year - 1)
  firstMod <- t0
  lastMod <- tn

  yrs <- sort(unique(datm$year[datm$rec == 1]))
  gaps <- if (length(yrs) > 1) diff(range(yrs)) else 1

  meta <- data.frame(species = species,
                     nRec_glob = ifelse(exists("nRec_glob"), nRec_glob, NA),
                     nRec_reg = ifelse(exists("nRec_reg"), nRec_reg, NA),
                     first = first, last = last,
                     gap = gaps,
                     firstMod = firstMod, lastMod = lastMod,
                     median = NA, P90 = NA,
                     visits_median = NA, visits_P90 = NA,
                     prop_list_one = NA, prop_repeats_grp = NA, prop_abs = NA,
                     EqualWt = NA, HighSpec = NA)

  return(list(raw_occ, meta))
}
