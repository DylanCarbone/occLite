#' @export
remove_bad_species <- function (Occ, threshold_sd, threshold_yrs, threshold_Rhat){
  if ("sd" %in% dimnames(Occ)[[3]] & "Rhat" %in% dimnames(Occ)[[3]]) {
    reliable <- Occ[, , "sd", drop = FALSE] < threshold_sd & 
      Occ[, , "Rhat", drop = FALSE] < threshold_Rhat
    reliable <- apply(reliable, c(1, 2), mean)
    reliable[reliable == 0] <- NA
    Occ <- Occ[, , !dimnames(Occ)[[3]] %in% c("sd", "Rhat"), 
               drop = FALSE]
  }
  else {
    sds <- apply(Occ, c(1, 2), sd, na.rm = TRUE)
    reliable <- sds < threshold_sd
    reliable[!reliable] <- NA
  }
  OccRel <- sapply(1:dim(Occ)[3], function(j) as.numeric(Occ[, , j]) * 
                     reliable, simplify = "array")
  OccRel <- OccRel[rowSums(reliable, na.rm = T) >= threshold_yrs, 
                   , , drop = FALSE]
  if (dim(OccRel)[1] == 0) 
    stop("None of your species meet the thresholds")
  attr(OccRel, "good_years") <- reliable[rowSums(reliable, 
                                                 na.rm = T) >= threshold_yrs, ]
  return(OccRel)
}
