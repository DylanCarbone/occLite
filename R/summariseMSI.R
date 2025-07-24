#' @title Summarise and Plot Occupancy Index Trends
#'
#' @description
#' A modified function to summarise and optionally visualise occupancy index results
#' from multi-species indicators (MSI). Produces either time-series plots of the 
#' index, bar plots of species trend assessments, or plots of index uncertainty 
#' versus number of contributing species. Optionally returns summary data tables instead of plots.
#'
#' @param label Character. A label to use as the plot title and file name if saving plots.
#' @param plotType Character. Type of plot to produce. Options are `"indicator"` (for time-series index), 
#' `"nSpecies"` (number of species through time), or `"uncertainty"` (width of interval vs number of species).
#' @param indicator A list containing the output of a multi-species indicator analysis, 
#' typically including components like `Summary`, `st`, `lt`, `final`, and `MetaData`.
#' @param method Character. The method used to compute trends, e.g. `"lambda"`. Required to enable trend bar plots.
#' @param plot Logical. If `TRUE`, plots are displayed; if `FALSE`, summary data is returned instead.
#' @param writePlot Logical. If `TRUE` and `plotType == "indicator"` with `method == "lambda"`, the plot will be saved to disk.
#' @param minYear Integer. The minimum year to include in plots.
#' @param maxYear Integer. The maximum year to include in plots.
#'
#' @return If `plot = TRUE`, returns a ggplot or grid of plots. If `plot = FALSE`, returns a list of:
#' \itemize{
#'   \item `out1`: A data.frame of summarised occupancy values by year
#'   \item `change`: A table of species trend categories across time periods
#'   \item `raw`: A data.frame of raw trend category labels for each species
#' }
#'
#' @examples
#' # Example using precomputed indicator
#' modified_summariseMSI(
#'   label = "Ants",
#'   plotType = "indicator",
#'   indicator = ants_indicator_result,
#'   method = "lambda",
#'   plot = TRUE,
#'   writePlot = FALSE,
#'   minYear = 1990,
#'   maxYear = 2020
#' )
#'
#' @export
modified_summariseMSI <- 
  function (label, plotType, indicator, method, plot = TRUE, writePlot = TRUE, minYear, maxYear)  
  {
    n <- max(indicator$Summary$Species_Number)
    nSpec <- indicator$Summary$Species_Number
    years <- minYear:maxYear
    if (plotType == "indicator") {
      p1 <- ggplot(data = indicator$Summary[indicator$Summary$year %in% years,], aes(x = years, y = indicator)) + 
        geom_ribbon(data = indicator$Summary[indicator$Summary$year %in% years,], aes(ymax = upper, 
                                                                                      ymin = lower), fill = "grey80") + 
        geom_line(size = 0.25) + geom_point(size = 0.25) + theme_linedraw() + ylab("Occupancy index") + 
        xlab("") + ylim(c(0, 150)) + ggtitle(label) + annotate("text", 
                                                               x = 1985, y = 30, label = paste(n, "species"))
      if (method == "lambda") {
        st <- indicator$st$species_assessment$category[!is.na(indicator$st$species_assessment$category)]
        st <- data.frame(st, rep(as.factor(1), length(st)))
        colnames(st) <- c("val", "type")
        lt <- indicator$lt$species_assessment$category[!is.na(indicator$lt$species_assessment$category)]
        lt <- data.frame(lt, rep(as.factor(2), length(lt)))
        colnames(lt) <- c("val", "type")
        dat <- rbind(st, lt)
        p2 <- ggplot(dat, aes(x = factor(type), fill = forcats::fct_rev(val))) + 
          geom_bar(position = "fill") + theme_linedraw() + 
          ylab("Proportion of species") + xlab("") + scale_x_discrete(labels = c("Short term", 
                                                                                 "Long term")) + guides(fill = guide_legend(title = ""))
        IndPlot <- grid.arrange(p1, p2, ncol = 2)
        print(IndPlot)
        if(writePlot){
          ggsave(IndPlot, filename = paste0("ind_", label, ".png"), height = 3, width = 9,
                 units = "in", dpi = 500)
        }
      }
      else {
        print(p1)
      }
    }
    else if (plotType == "nSpecies") {
      p1 <- ggplot(data = NULL, aes(x = years, y = indicator$Summary$Species_Number)) + 
        geom_line() + geom_point() + theme_linedraw() + ylab("Number of species contributing") + 
        xlab("") + ggtitle(label)
      print(p1)
    }
    else {
      width <- summary$upper - summary$lower
      p1 <- ggplot(data = NULL, aes(y = width, x = ind$summary$Species_Number)) + 
        geom_point() + theme_linedraw() + ylab("Width of credible interval") + 
        xlab("Number of species contributing") + ggtitle(label)
      print(p1)
    }
    if (!plot) {
      indicator$Summary$year <- indicator$Summary$year + (minYear - 
                                                            1)
      out1 <- indicator$Summary
      changeLT <- data.frame(table(indicator$MetaData$species_change$category))
      rawLT <- indicator$MetaData$species_change[, 1]
      L <- length(rawLT)
      changeST <- data.frame(table(indicator$st$species_assessment$category))[, 
                                                                              2]
      rawST <- indicator$st$species_assessment[, 1]
      rawST <- c(rawST, rep(NA, (L - length(rawST))))
      changeFinal <- data.frame(table(indicator$final$species_assessment$category))[, 
                                                                                    2]
      rawFinal <- indicator$final$species_assessment[, 1]
      rawFinal <- c(rawFinal, rep(NA, (L - length(rawFinal))))
      change <- data.frame(changeLT, changeST, changeFinal)
      colnames(change) <- c("Trend", "Long term", "Short term", 
                            "Final year")
      raw <- data.frame(rawLT, rawST, rawFinal)
      colnames(raw) <- c("Long term", "Short term", "Final year")
      out <- list(out1, change, raw)
    }
    else {
      out <- p1
    }
  }