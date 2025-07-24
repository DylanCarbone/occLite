#' \code{loadRFile} - function to load Rdata
#' 
#' @description This function loads sparta Rdata object
#'
#' @param fileName name of the file to load
#' 
#' @return file loaded

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  return(get(ls()[ls() != "fileName"]))
}
