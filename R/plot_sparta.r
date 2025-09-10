#' Plot occDet Objects
#' 
#' @param occDet An object of class occDet
#' @param main The plot title, defaults to the species name
#' @param reg_agg The name of a region or region aggregate to plot.
#' If '' (default) then the overall occupancy estimates are plotted
#' 
#' @import ggplot2
#' @export
plot_sparta <- function(occDet, main = occDet$SPP_NAME, reg_agg = NULL) {
  
  # gets summary output from the BUGS files 
  spp_data <- as.data.frame(occDet$BUGSoutput$summary)
  
  if (!is.null(reg_agg)) {
    reg_agg <- paste0(".r_", reg_agg)
  } else {
    reg_agg <- ""
  }
  
  # take psi.fs rows = yearly proportion of occupied cells
  spp_data$occDet <- row.names(spp_data)
  new_data <- spp_data[grepl(paste0("^psi.fs", reg_agg, "\\["), spp_data$occDet), ]
  new_data$year <- (occDet$min_year - 1) + 
    as.numeric(gsub(paste0("psi.fs", reg_agg), "", gsub("\\[|\\]","", row.names(new_data))))
  
  # rename columns
  names(new_data) <- gsub("2.5%", "quant_025", names(new_data))
  names(new_data) <- gsub("97.5%", "quant_975", names(new_data))
  
  # add Rhat quality flag
  new_data$rhat_threshold <- ifelse(new_data$Rhat < 1.1, "Good (<1.1)", "Bad (>1.1)")
  
  # plot (occupancy scale 0–1)
  p = ggplot(new_data, aes(x = year, y = mean)) +
    theme_minimal() +
    geom_ribbon(aes(ymin = quant_025, ymax = quant_975), fill = "red", alpha = 0.15) +
    geom_line(colour = "darkred", size = 1.2) +
    geom_point(aes(colour = rhat_threshold), size = 3) +
    xlab("Year") +
    ylab("Occupancy") +
    ggtitle(main) +
    ylim(0, 1) +
    theme(
      plot.title = element_text(lineheight = .8, face = "bold"),
      legend.position = "bottom"
    ) +
    scale_colour_manual(
      name = "Rhat",
      values = c("Good (<1.1)" = "black", "Bad (>1.1)" = "orange", "start" = "blue", "end" = "red")
    )

  return(p)
}
