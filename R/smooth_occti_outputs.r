#' Generate OCCTI Occupancy Plots
#'
#' Generates LOESS-smoothed occupancy plots for each species and optionally saves them to disk.
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
smooth_occti_outputs <- function(occti_outputs,
                                 save_dir = "occupancy_plots",
                                 save_plots = TRUE,
                                 span = 0.75,
                                 plot_width = 18,
                                 plot_height = 6,
                                 n_iter = 1000) {
                             
  
  # Create output folder if saving is enabled and directory does not exist
  if (save_plots && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Prepare list to hold LOESS summaries for each species and occupancy plots
  loess_summaries <- list()
  occupancy_plots <- list()
  
  # Loop through each species in the input list
  for (species in names(occti_outputs)) {
    
    # Extract the occupancy index data, bound to be less than 1 and more than 0, then logit-transform
    occ_data <- occti_outputs[[species]]$Index %>%
      mutate(psiA = logit(bound_numbers(psiA)),
      psiA_L = logit(bound_numbers(psiA_L)), 
      psiA_U = logit(bound_numbers(psiA_U)))
    
    # Generate simulated occupancy values using normal distribution
    psiA_draws <- do.call(rbind, lapply(1:nrow(occ_data), function(i) {
      year <- occ_data$Year[i]
      mean_val <- occ_data$psiA[i]
      # Estimate standard deviation from 95% CI bounds
      sd_val <- (occ_data$psiA_U[i] - occ_data$psiA_L[i]) / (2 * 1.96)
      
      # Simulate n_iter draws for the year
      data.frame(
        Year = year,
        Iteration = 1:n_iter,
        Simulated_psiA = rnorm(n_iter, mean = mean_val, sd = sd_val)
      )
    }))
    
    # Fit LOESS model per iteration and predict occupancy per year
    loess_predictions <- psiA_draws %>%
      group_by(Iteration) %>%
      do({
        # Try to fit a LOESS curve to each simulated iteration
        mod <- try(loess(Simulated_psiA ~ Year, data = ., span = span), silent = TRUE)
        pred_vals <- if (inherits(mod, "try-error")) rep(NA, nrow(.)) else predict(mod, newdata = data.frame(Year = .$Year))
        data.frame(Year = .$Year, Pred = pred_vals)
      }) %>%
      ungroup()
    
    # Summarise LOESS results across all iterations
    loess_summary <- loess_predictions %>%
      group_by(Year) %>%
      summarise(
        psiA_loess_mean = mean(Pred, na.rm = TRUE),
        psiA_loess_lower = quantile(Pred, 0.025, na.rm = TRUE),
        psiA_loess_upper = quantile(Pred, 0.975, na.rm = TRUE),
        psiA_loess_se = sd(Pred, na.rm = TRUE) / sqrt(sum(!is.na(Pred))),
        .groups = "drop"
      )
    
    # Store the summary in the output list
    loess_summaries[[species]] <- loess_summary

    # Create the final ggplot with original and smoothed estimates
    p <- ggplot(occ_data, aes(x = Year)) + 
      geom_line(data = loess_summary, aes(y = inv.logit(psiA_loess_mean)), colour = "darkred", size = 1.2) +
      geom_ribbon(data = loess_summary, aes(ymin = inv.logit(psiA_loess_lower), ymax = inv.logit(psiA_loess_upper)), fill = "red", alpha = 0.15) +
      labs(x = "Year", y = "Occupancy Index") +
      theme_minimal() +
      ggtitle(paste(species, "- loess span =", span)) +
      ylim(0, 1)
    
    occupancy_plots[[species]] = p
    
    # Save plot if enabled
    if (save_plots) {
      file_name <- file.path(save_dir, paste0(species, ".png"))
      ggsave(filename = file_name, plot = p, width = plot_width, height = plot_height, dpi = 300)
    }
  }
  
  # Return the list of LOESS summaries and the occupancy plots
  return(list("loess_summaries" = loess_summaries, "occupancy_plots" = occupancy_plots))
}
