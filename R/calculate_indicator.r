#' Calculate biodiversity indicator using lambda or BMA method
#'
#' This function calculates a biodiversity indicator from occti predictions
#' using either the lambda or BMA (Bayesian meta analysis) approach.
#' It filters data by year, transforms predictions as needed, and returns the full
#' indicator object, a summary table, and a final trend assessment.
#'
#' @param loess_predictions A data frame containing LOESS model predictions with columns:
#'   \code{species}, \code{year}, and \code{pred}.
#' @param method A character string indicating the indicator calculation method. Must be one of:
#'   \code{"lambda"} or \code{"bma"}.
#' @param min_year Optional numeric value specifying the minimum year to include. If \code{NULL},
#'   the minimum year in \code{loess_predictions} will be used.
#' @param max_year Optional numeric value specifying the maximum year to include. If \code{NULL},
#'   the maximum year in \code{loess_predictions} will be used.
#' @param bma_ind Character string for selecting the BMA indicator type, only used when \code{method = "bma"}.
#'   Defaults to \code{NULL}.
#'
#' @return A list with three elements:
#'   \item{indicator}{The full output from either \code{lambda_indicator()} or \code{BRCindicators::bma()}.}
#'   \item{summary}{A data frame with yearly indicator values, confidence intervals, species count, and year.}
#'   \item{final}{A data frame summarising trend assessment statistics based on the indicator.}
#'
#' @details
#' For method \code{"bma"}, predictions are transformed to the logit scale and averaged across iterations
#' for each species-year combination. Standard errors are calculated and passed to
#' \code{BRCindicators::bma()}.
#'
#' For method \code{"lambda"}, predictions are reshaped into a 3D array using \code{reshape2::acast()},
#' and passed to \code{lambda_indicator()}. No transformation is applied, but quantiles are calculated
#' from posterior distributions across iterations.
#'
#' @seealso \code{\link{lambda_indicator}}, \code{BRCindicators::bma}, \code{\link{trend_assessment}}
#'
#' @import dplyr BRCindicators
#' @importFrom reshape2 acast
#' @export

calculate_indicator <- function(loess_predictions, sparta_output = TRUE, method, min_year = NULL, max_year = NULL, bma_ind = NULL){
# Pivot the sparta outputs (if the outputs are sparta)
if(sparta_output){
   loess_predictions = translate_sparta(loess_predictions)
}

if (!method %in% c("lambda", "bma")){ 
    stop("Method must be one of lambda or bma")}

min_year <- ifelse(is.null(min_year), 
            yes = min(loess_predictions$year),
            no = startYear)
max_year <- ifelse(is.null(max_year), 
            yes = max(loess_predictions$year),
            no = endYear)

loess_predictions = loess_predictions %>% filter(year >= min_year, year <= max_year)

if(method == "bma"){
    
mean_se_logit <- loess_predictions %>%
    mutate(pred = logit(bound_for_logit(pred))) %>%
    group_by(species, year) %>%
    summarise(index = mean(pred),
    se = sd(pred, na.rm = TRUE) / sqrt(sum(!is.na(pred))))

indicator_output <- BRCindicators::bma(data = mean_se_logit, seFromData = TRUE, 
                            m.scale = "logit") # keep at loge for now

if (bma_ind != "prime") {
    summary <- data.frame(indicator = indicator_output$Index.M, 
                        lower = indicator_output$lowerCI.M, 
                        upper = indicator_output$upperCI.M, 
                        Species_Number = length(unique(mean_se_logit$species)),
                        year = indicator_output$Year)
} else {
    summary <- data.frame(indicator = indicator_output$Index.Mprime, 
                        lower = indicator_output$lowerCI.Mprime, 
                        upper = indicator_output$upperCI.Mprime, 
                        Species_Number = length(unique(mean_se_logit$species)), # should this be the number of species with data each year?
                        year = indicator_output$Year)
}

trend_assessment_summary <- trend_assessment(dat = indicator_output, summary = summary, method = method, 
                                    start_year = min(indicator_output$Year), end_year = max(indicator_output$Year))
}

else{

# Create the array using acast
occ_array <- reshape2::acast(loess_predictions, species ~ year ~ iteration, value.var = "pred")

# Calculate lambda
indicator_output <- lambda_indicator(input = occ_array, index = 100, 
                                    threshold_sd = Inf, threshold_Rhat = Inf, threshold_yrs = 10, 
                                    upperQuantile = 0.975, lowerQuantile = 0.025)

summary <- indicator_output$summary

year = indicator_output$summary %>%
  filter(!is.na(indicator)) %>%
  summarise(max_year = max(year)) %>%
  pull(max_year)

trend_assessment_summary <- trend_assessment(
    indicator_output, summary = summary, method = method, start_year = year, end_year = year)  
}

return(list("indicator" = indicator_output, "summary" = summary, "final" = trend_assessment_summary))

}