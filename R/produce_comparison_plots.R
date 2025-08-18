#' Function to produce comparison plots for DEFRA report
#' @import ggplot2 patchwork
#' @export
#' 
produce_comparison_plots <- function(occti_loess_samples, sparta_occupancy_estimates_folder, output_folder = "sparta_occti_comparison", sparta_region){

    sparta_occ_species = gsub(".rdata", "", basename(list.files(sparta_occupancy_estimates_folder)))
    
    for (species_occ in unique(occti_loess_samples$species)){

        if (!species_occ %in% sparta_occ_species){
            warning(paste("species", species_occ, "is missing from sparta outputs. Skipping..."))
            next
        }

        loess_summary_species <- occti_loess_samples %>%
        filter(species == species_occ) %>%
        group_by(year) %>%
        summarise(
            psiA_loess_mean = mean(pred, na.rm = TRUE),
            psiA_loess_lower = quantile(pred, 0.025, na.rm = TRUE),
            psiA_loess_upper = quantile(pred, 0.975, na.rm = TRUE)
        )

        # Plot
        plot_occti <- ggplot(loess_summary_species, aes(x = year)) + 
        geom_line(data = loess_summary_species, aes(y = psiA_loess_mean), colour = "darkred", size = 1.2) +
        geom_ribbon(data = loess_summary_species, aes(ymin = psiA_loess_lower, ymax = psiA_loess_upper), fill = "red", alpha = 0.15) +
        labs(x = "Year", y = "Occupancy") +
        theme_minimal() +
        ylim(0,1)

        # Obtain sparta outputs
        load(file.path(sparta_occupancy_estimates_folder, paste0(species_occ, ".rdata")))

        plot_sparta <- plot_sparta(out, reg_agg = sparta_region, main = NULL)

        combined = plot_sparta + plot_occti +
        plot_annotation(
            title = paste(species_occ, "sparta and occti comparison")
        )

        if(!file.exists(output_folder)){
            dir.create(output_folder)
        }

        # Save to file
        ggsave(file.path(output_folder, paste0(species_occ, "_comparison.png")), plot = combined, width = 12, height = 4, dpi = 300)
    }

    return(NULL)
}