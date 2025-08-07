#' @export
trend_assessment <- function (dat, summary, method = "lambda", start_year, end_year, 
                                       species_stat = "mean"){

  if (!method %in% c("lambda", "bma")) 
    stop("Method must be one of 'lambda' or 'bma'")
  if (method == "lambda") {

    sp_assess <- BRCindicators:::species_assessment(dat = dat$LogLambda, 
                                                    method = method, start_year = start_year, end_year = end_year, # LB replace start_year = start_year +1 with start_year =start_year
                                                    species_stat = species_stat, plot = FALSE)
    ind_assessment <- BRCindicators:::indicator_assessment(summary_table = summary, 
                                                           start_year = start_year, end_year = end_year)
    return(list(species_assessment = sp_assess, indicator_asssessment = ind_assessment)) 
  }
  else {
    # the column name has to be year (lowercase) inside species_assessment
    colnames(dat) <-gsub("Year", "year", colnames(dat))
    
    sp_assess <- BRCindicators:::species_assessment(dat = dat, method = method, 
                                                    start_year = start_year, end_year = end_year, species_stat = species_stat, 
                                                    plot = FALSE)
    ind_assessment <- BRCindicators:::indicator_assessment(summary_table = summary, 
                                                           start_year = start_year, end_year = end_year)
    return(list(species_assessment = sp_assess, indicator_asssessment = ind_assessment))
  }
}