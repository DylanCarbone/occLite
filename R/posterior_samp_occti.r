#' Generate OCCTI posterior samples
#'
#' Generates smoothed posterior samples of occupancy from loess model predictions. Optionally saves models to disk.
#' Also returns a list of LOESS-smoothed values (mean and bounds) per species.
#'
#' @param occti_outputs A named list of occupancy outputs as returned by `load_occti_outputs()`.
#' @param save_dir A character string for the directory where plots should be saved (if `save_plots = TRUE`).
#' @param save_plots Logical. If TRUE, plots will be saved to `save_dir`.
#' @param span Numeric. The smoothing parameter passed to the `loess()` function (e.g. 0.75).
#' @param plot_width Numeric. Width of saved plots in inches. Used only if `save_plots = TRUE`.
#' @param plot_height Numeric. Height of saved plots in inches. Used only if `save_plots = TRUE`.
#' @param n_iter Integer. Number of iterations to use to sample within mean and standard deviation`.
#'
#' @return A named list where each element is a data frame of LOESS-smoothed occupancy values for one species.
#'
#' @import ggplot2 dplyr patchwork boot
#' @export
posterior_samp_occti <- function(occti_outputs,
                                 save_dir = "occupancy_plots",
                                 save_plots = TRUE,
                                 span = 0.75,
                                 plot_width = 18,
                                 plot_height = 6,
                                 n_iter = 1000) {
  
  if (save_plots && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  loess_predictions_all <- list()
  
  for (species in names(occti_outputs)) {
    
    occ_data <- occti_outputs[[species]]$Index
    
    psiA_draws <- do.call(rbind, lapply(1:nrow(occ_data), function(i) {
      year <- occ_data$Year[i]
      mean_val <- occ_data$psiA[i]
      sd_val <- (occ_data$psiA_U[i] - occ_data$psiA_L[i]) / (2 * 1.96)
      
      data.frame(
        year = year,
        iteration = 1:n_iter,
        simulated_psiA = rnorm(n_iter, mean = mean_val, sd = sd_val)
      )
    })) %>%
      filter(simulated_psiA <= 1, simulated_psiA >= 0)
    
    years_range <- range(occ_data$Year)
    years <- years_range[1]:years_range[2]
    
    loess_predictions <- suppressMessages(suppressWarnings(psiA_draws %>%
      group_by(iteration) %>%
      do({
        mod <- try(loess(simulated_psiA ~ year, data = ., span = span), silent = TRUE)
        if (inherits(mod, "try-error")) {
          data.frame(year = years, pred = NA)
        } else {
          pred <- try(predict(mod, newdata = data.frame(year = years)), silent = TRUE)
          if (inherits(pred, "try-error")) {
            data.frame(year = years, pred = NA)
          } else {
            data.frame(year = years, pred = pred)
          }
        }
      }) %>%
      ungroup() %>%
      filter(!is.na(pred)) %>%
      mutate(pred = bound_zero_one(pred), species = species)))
        
    # Only store if we have valid iterations for this species
    if (nrow(loess_predictions) > 0) {
      loess_predictions_all[[species]] <- loess_predictions
      
      if (save_plots) {
        loess_summary_occ <- loess_predictions %>%
          group_by(year) %>%
          summarise(
            psiA_loess_mean = mean(pred, na.rm = TRUE),
            psiA_loess_lower = quantile(pred, 0.025, na.rm = TRUE),
            psiA_loess_upper = quantile(pred, 0.975, na.rm = TRUE),
            psiA_loess_se = sd(pred, na.rm = TRUE) / sqrt(sum(!is.na(pred))),
            .groups = "drop"
          )
        
        p <- ggplot(occ_data, aes(x = year)) +
          geom_line(data = loess_summary_occ, aes(y = psiA_loess_mean), colour = "darkred", size = 1.2) +
          geom_ribbon(data = loess_summary_occ, aes(ymin = psiA_loess_lower, ymax = psiA_loess_upper),
                      fill = "red", alpha = 0.15) +
          labs(x = "Year", y = "Occupancy Index") +
          theme_minimal() +
          ggtitle(paste(species, "- loess span =", span)) +
          ylim(0, 1)
        
        file_name <- file.path(save_dir, paste0(species, ".png"))
        ggsave(filename = file_name, plot = p, width = plot_width, height = plot_height, dpi = 300)
      }
    } else {
      stop(sprintf("No valid iterations for species '%s'. This was likely because there were insufficient years to fit the loess model. Please filter the results first using occLite::filter_occupancy()", species))
    }
  }
  
  # Remove NULLs before binding
  loess_predictions_all <- loess_predictions_all[!sapply(loess_predictions_all, is.null)]
  
  # Combine into a single dataframe
  loess_predictions_df <- bind_rows(loess_predictions_all)
  
  return(loess_predictions_df)
}
