#' @export
lambda_indicator <- function (input, index = 100, threshold_sd = 0.2, threshold_Rhat = 1.1, 
                                       threshold_yrs = 20, upperQuantile = 0.975, lowerQuantile = 0.025, 
                                       sample_size = NULL, year_range = NULL, region = NULL){
  Occ <- BRCindicators:::getData(input = input, sample_size = sample_size, 
                                 region = region)
  if (!is.null(year_range)){ 
    Occ <- subset_years(Occ, year_range)}
  nsp1 <- dim(Occ)[1]
  Occ <- remove_bad_species(Occ = Occ, threshold_sd = threshold_sd, 
                                     threshold_yrs = threshold_yrs, threshold_Rhat = threshold_Rhat)
  good_years <- attr(Occ, "good_years")
  nsp2 <- dim(Occ)[1]
  if (dim(Occ)[1] == 0) {
    stop("No species meet the threshold criteria")
  }
  else if ((nsp1 - nsp2) > 0) {
    message(paste(nsp1 - nsp2, "Species have been removed as they don't meet your thresholds"))
  }
  Occ <- car::logit(Occ, adjust = 0.001)
  LogLambda <- apply(Occ, c(1, 3), BRCindicators:::lambda_calc)
  LogLambda <- aperm(LogLambda, c(2, 1, 3))
  for (i in 1:nrow(LogLambda)) {
    firstZero <- match(0, LogLambda[i, , ])
    LogLambda[i, firstZero, ] <- NA
  }
  dimnames(LogLambda) <- dimnames(Occ)
  Delta <- apply(LogLambda, c(2, 3), mean, na.rm = T)
  Delta[1, ] <- 0
  Theta <- matrix(data = NA, nrow = nrow(Delta), ncol = ncol(Delta), 
                  dimnames = dimnames(Delta))
  Theta[!is.na(rowMeans(Delta)), ] <- apply(Delta[!is.na(rowMeans(Delta)), 
  ], 2, cumsum)
  Lambda <- exp(Theta)
  Indicator_data <- index * Lambda
  indicator <- apply(Indicator_data, 1, mean)
  Indicator_CI <- t(apply(Indicator_data, 1, quantile, probs = c(lowerQuantile, 
                                                                 upperQuantile), na.rm = TRUE))
  colnames(Indicator_CI) <- c("lower", "upper")
  summary_table <- as.data.frame(cbind(indicator, Indicator_CI))
  summary_table$year <- as.numeric(row.names(summary_table))
  if (any(is.na(summary_table))) 
    warning("Data not available for all years in output due to threshold data removal")
  sp_change <- BRCindicators:::species_assessment(LogLambda, start_year = min(as.numeric(dimnames(LogLambda)[[2]])) + 
                                                    1, end_year = max(as.numeric(dimnames(LogLambda)[[2]])))
  LogLambdaPart <- LogLambda[, , 1]
  LogLambdaPartTF <- (!(is.na(LogLambdaPart)))
  NperYear <- colSums(LogLambdaPartTF)
  summary_table$Species_Number <- NperYear
  return(list(summary = summary_table, LogLambda = LogLambda, 
              ind_data = Indicator_data, species_change = sp_change, 
              good_years = good_years))
}
