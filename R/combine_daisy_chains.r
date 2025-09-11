#' Combine Daisy-Chained Model Outputs Across Chains
#'
#' For each species, loads only the *final* iteration chunk across
#' all chains, combines them, and trims to the last 2000 posterior draws.
#' Metadata (regions, region_aggs, etc.) is retrieved from the *first* chunk.
#'
#' @param input_path Character. Directory with .rdata chain files.
#' @param output_path Character. Directory to save combined outputs (default = "combined_chains").
#' @param keep_n Integer. Number of posterior draws to keep from the end.
#'   Default = 2000.
#'
#' @return Invisibly returns vector of saved file paths.
#'
#' @import dplyr
#' @import stringr
#'
#' @export
#'
combine_daisy_chains <- function(input_path,
                                 output_path = "combined_chains",
                                 keep_n = 2000) {

  files <- list.files(input_path, pattern = "\\.rdata$", full.names = TRUE)
  if (length(files) == 0) stop("No .rdata files found in ", input_path)

  # Parse species, iteration, chain from file names
  df <- tibble(file = files) %>%
    mutate(base = basename(file) %>% str_remove("\\.rdata$")) %>%
    tidyr::separate(base,
                    into = c("sp_group", "sp_id", "iteration", "chain"),
                    sep = "_", convert = TRUE) %>%
    mutate(species = paste(sp_group, sp_id, sep = "_"))

  outfiles <- c()

  for (sp in unique(df$species)) {
    message("Combining species: ", sp)

    df_sp <- df %>% filter(species == sp)

    # --- Load metadata from the *first* chunk ---
    min_iter <- min(df_sp$iteration, na.rm = TRUE)
    first_file <- df_sp %>% filter(iteration == min_iter) %>% slice(1) %>% pull(file)
    load(first_file)  # loads object "out"
    out_template <- out   # keep metadata from first chunk

    # --- Merge chains from the *final* chunk ---
    max_iter <- max(df_sp$iteration, na.rm = TRUE)
    df_final <- df_sp %>% filter(iteration == max_iter) %>% arrange(chain)

    sims_list <- NULL
    for (f in df_final$file) {
      load(f)  # loads "out"
      sims_chain <- out$BUGSoutput$sims.list
      sims_list <- if (is.null(sims_list)) sims_chain else Map(rbind, sims_list, sims_chain)
    }

    # --- Take the last keep_n draws across all chains ---
    sims_list <- lapply(sims_list, function(mat) {
      n <- nrow(mat)
      start <- max(1, n - keep_n + 1)
      mat[start:n, , drop = FALSE]
    })

    n.sims_total <- nrow(sims_list[[1]])

    # --- Replace BUGSoutput sims with trimmed version ---
    out_template$BUGSoutput$sims.list <- sims_list
    out_template$BUGSoutput$n.sims   <- n.sims_total

    message(paste("kept", n.sims_total, "posterior draws for", sp))

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
