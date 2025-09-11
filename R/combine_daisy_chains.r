#' Combine Daisy-Chained Model Outputs Across Iteration Chunks and Chains
#'
#' Combines .rdata files from daisy-chained model runs (multiple iteration chunks,
#' multiple chains per chunk) into a single model object per species.
#'
#' @param input_path Character. Directory with .rdata chain files.
#' @param output_path Character. Directory to save combined outputs (default = "combined_chains").
#'
#' @return Invisibly returns vector of saved file paths.
#' 
#' @import dplyr
#' @import stringr
#'
#' @export
#'
combine_daisy_chains <- function(input_path, output_path = "combined_chains") {

  files <- list.files(input_path, pattern = "\\.rdata$", full.names = TRUE)
  if (length(files) == 0) stop("No .rdata files found in ", input_path)

  # Parse species, iteration, chain from file names
  df <- tibble(file = files) %>%
    mutate(base = basename(file) %>% str_remove("\\.rdata$")) %>%
    tidyr::separate(base, into = c("sp_group", "sp_id", "iteration", "chain"),
                    sep = "_", convert = TRUE) %>%
    mutate(species = paste(sp_group, sp_id, sep = "_"))

  outfiles <- c()

  for (sp in unique(df$species)) {
    message("Combining species: ", sp)

    df_sp <- df %>%
      filter(species == sp) %>%
      arrange(iteration, chain)

    sims_list <- NULL
    n.sims_total <- 0
    out_template <- NULL

    for (iter_val in unique(df_sp$iteration)) {
      chunk_files <- df_sp %>% filter(iteration == iter_val)

      sims_chunk <- NULL
      for (f in chunk_files$file) {
        load(f)  # loads object "out"
        if (is.null(out_template)) out_template <- out
        sims_chain <- out$BUGSoutput$sims.list

        sims_chunk <- if (is.null(sims_chunk)) sims_chain else Map(rbind, sims_chunk, sims_chain)

        # increment per file
        n.sims_total <- n.sims_total + out$BUGSoutput$n.sims
      }

      sims_list <- if (is.null(sims_list)) sims_chunk else Map(rbind, sims_list, sims_chunk)
    }

    # Replace sims in template
    out_template$BUGSoutput$sims.list <- sims_list
    out_template$BUGSoutput$n.sims   <- n.sims_total

    message(paste("for species", sp, "there were", n.sims_total, "total posterior draws"))

    # Ensure output dir exists
    if (!dir.exists(output_path)) {
      dir.create(output_path, recursive = TRUE)
    }

    outfile <- file.path(output_path, paste0(sp, ".rdata"))
    out <- out_template
    save(out, file = outfile)

    outfiles <- c(outfiles, outfile)
  }

  invisible(outfiles)
}
