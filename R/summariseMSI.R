#' Summarise and Plot Multi-Species Occupancy Indicators
#'
#' Creates summary plots for a multi-species occupancy indicator, including:
#' \itemize{
#'   \item A time-series plot of the occupancy index with uncertainty intervals.
#'   \item A stacked bar plot showing the proportion of species in each final-year trend category.
#' }
#' Optionally saves the combined plot to disk.
#'
#' @param indicator List. Output from a multi-species indicator analysis in 
#'   \code{sparta_indicator} format. Must contain:
#'   \itemize{
#'     \item \code{indicator$summary} — Data frame with yearly index values, 
#'           uncertainty bounds, and number of contributing species.
#'     \item \code{indicator$final$species_assessment} — Data frame with 
#'           per-species trend assessment categories (factor or character).
#'   }
#' @param label Character. A label for the plots; used as the plot title and 
#'   in the saved file name (if \code{writePlot = TRUE}).
#' @param writePlot Logical. If \code{TRUE}, saves the combined plot to a 
#'   PNG file in the working directory with filename \code{"ind_<label>.png"}.
#' @param minYear Integer. First year to include in the time-series plot.
#' @param maxYear Integer. Last year to include in the time-series plot.
#'
#' @return A named list with two ggplot objects:
#' \describe{
#'   \item{\code{indicator}}{Time-series occupancy index plot with shaded uncertainty ribbon.}
#'   \item{\code{species_change}}{Final-year stacked bar plot showing the proportion of species
#'     in each trend category.}
#' }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Extracts the time-series summary data from \code{indicator$summary} and filters to the 
#'         specified year range (\code{minYear}–\code{maxYear}).
#'   \item Plots the occupancy index with a ribbon showing the 95\% credible interval, 
#'         points and lines for yearly estimates, and the number of species annotated on the plot.
#'   \item Extracts the final-year trend categories from 
#'         \code{indicator$final$species_assessment$category} and plots them as proportions.
#'   \item Combines the two plots side-by-side using \pkg{patchwork}.
#'   \item Saves the combined plot to file if \code{writePlot = TRUE}.
#' }
#'
#' The y-axis limit for the time-series plot is set slightly above the maximum upper bound 
#' to provide visual padding. Trend category bars are ordered in reverse to match the intended legend order.
#'
#' @import patchwork
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point theme_linedraw ylab xlab ylim ggtitle annotate geom_bar guides guide_legend ggsave
#' @importFrom forcats fct_rev
#'
#' @export

summariseMSI <- function(indicator, label, writePlot = TRUE, minYear, maxYear)  {

    summary_df <- indicator$summary
    n <- max(summary_df$Species_Number)
    years <- minYear:maxYear

    plot_df = summary_df[summary_df$year %in% years,]
    p1 <- ggplot(plot_df, aes(x = year, y = indicator)) + 
      geom_ribbon(aes(ymax = upper, ymin = lower), fill = "grey80") + 
      geom_line(size = 0.25) + 
      geom_point(size = 0.25) + 
      theme_linedraw() + 
      ylab("Occupancy index") + 
      xlab("") + 
      ylim(c(0, max(plot_df$upper + 50))) + 
      ggtitle(label) + 
      annotate("text", x = min(years) + 1, y = 30, label = paste(n, "species"))

    final <- indicator$final$species_assessment$category
    final <- final[!is.na(final)]
    final_df <- data.frame(val = final, type = factor("Final year", levels = "Final year"))

    p2 <- ggplot(final_df, aes(x = type, fill = forcats::fct_rev(val))) + 
      geom_bar(position = "fill") + 
      theme_linedraw() + 
      ylab("Proportion of species") + 
      xlab("") + 
      guides(fill = guide_legend(title = ""))

    output <- list(indicator = p1, species_change = p2)
    combined_plot <- p1 + p2

    if (writePlot) {
      ggsave(combined_plot, filename = paste0("ind_", label, ".png"), height = 3, width = 9, units = "in", dpi = 500)
    }

    print(combined_plot)

    return(output)

  }