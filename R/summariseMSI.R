#' @title Summarise and Plot Occupancy Index Trends
#'
#' @description
#' A modified function to summarise and optionally visualise occupancy index results
#' from multi-species indicators (MSI). Produces either time-series plots of the 
#' index, bar plots of species trend assessments, or plots of index uncertainty 
#' versus number of contributing species. Optionally returns summary data tables instead of plots.
#'
#' This version is adapted for `sparta_indicator`-style output, where the main summary
#' data are stored in `indicator$indicator$summary`, and trend assessments are found
#' in `indicator$final$species_assessment`.
#'
#' @param label Character. A label to use as the plot title and file name if saving plots.
#' @param plotType Character. Type of plot to produce. Options are "indicator" (for time-series index), 
#' "nSpecies" (number of species through time), or "uncertainty" (width of interval vs number of species).
#' @param indicator A list containing the output of a multi-species indicator analysis. Must include
#'   elements named `indicator$indicator$summary` and `indicator$final$species_assessment`.
#' @param method Character. The method used to compute trends, e.g. "lambda". Required to enable trend bar plots.
#' @param plot Logical. If TRUE, plots are displayed; if FALSE, summary data is returned instead.
#' @param writePlot Logical. If TRUE and plotType == "indicator" with method == "lambda", the plot will be saved to disk.
#' @param minYear Integer. The minimum year to include in plots.
#' @param maxYear Integer. The maximum year to include in plots.
#'
#' @return If plot = TRUE, returns a ggplot or grid of plots. If plot = FALSE, returns a list of:
#' \itemize{
#'   \item out1: A data.frame of summarised occupancy values by year
#'   \item change: A table of species trend categories across time periods
#'   \item raw: A data.frame of raw trend category labels for each species
#' }
#'
#' @export
#' 
summariseMSI <- function(label, plotType, indicator, method, plot = TRUE, writePlot = TRUE, minYear, maxYear)  {
    summary_df <- indicator$indicator$summary
    n <- max(summary_df$Species_Number)
    years <- minYear:maxYear

    if (plotType == "indicator") {
      p1 <- ggplot(summary_df[summary_df$year %in% years,], aes(x = year, y = indicator)) + 
        geom_ribbon(aes(ymax = upper, ymin = lower), fill = "grey80") + 
        geom_line(size = 0.25) + 
        geom_point(size = 0.25) + 
        theme_linedraw() + 
        ylab("Occupancy index") + 
        xlab("") + 
        ylim(c(0, 150)) + 
        ggtitle(label) + 
        annotate("text", x = min(years) + 1, y = 30, label = paste(n, "species"))

      if (method == "lambda") {
        final <- indicator$final$species_assessment$category
        final <- final[!is.na(final)]
        final_df <- data.frame(val = final, type = factor("Final year", levels = "Final year"))

        p2 <- ggplot(final_df, aes(x = type, fill = forcats::fct_rev(val))) + 
          geom_bar(position = "fill") + 
          theme_linedraw() + 
          ylab("Proportion of species") + 
          xlab("") + 
          guides(fill = guide_legend(title = ""))

        IndPlot <- gridExtra::grid.arrange(p1, p2, ncol = 2)
        print(IndPlot)

        if (writePlot) {
          ggsave(IndPlot, filename = paste0("ind_", label, ".png"), height = 3, width = 9, units = "in", dpi = 500)
        }
      } else {
        print(p1)
      }
    } else if (plotType == "nSpecies") {
      p1 <- ggplot(summary_df[summary_df$year %in% years,], aes(x = year, y = Species_Number)) + 
        geom_line() + 
        geom_point() + 
        theme_linedraw() + 
        ylab("Number of species contributing") + 
        xlab("") + 
        ggtitle(label)
      print(p1)
    } else if (plotType == "uncertainty") {
      width <- summary_df$upper - summary_df$lower
      p1 <- ggplot(summary_df, aes(x = Species_Number, y = width)) + 
        geom_point() + 
        theme_linedraw() + 
        ylab("Width of credible interval") + 
        xlab("Number of species contributing") + 
        ggtitle(label)
      print(p1)
    }

    if (!plot) {
      summary_df$year <- summary_df$year + (minYear - 1)
      out1 <- summary_df
      rawFinal <- indicator$final$species_assessment[, 1]
      changeFinal <- data.frame(table(indicator$final$species_assessment$category))[, 2]
      change <- data.frame(Final = changeFinal)
      raw <- data.frame(Final = rawFinal)
      out <- list(out1, change, raw)
    } else {
      out <- p1
    }

    return(out)
  }
