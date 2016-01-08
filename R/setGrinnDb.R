#'Set Grinn database location
#'@description set Grinn database location for the currently working environment.
#'@usage setGrinnDb()
#'@param url of Grinn database location
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@export
#'@examples
#'# Set Grinn database location to local server
#'setGrinnDb("http://localhost:7474/db/data/")
setGrinnDb <- function(url){
  assign("nld", url, envir = .GlobalEnv)
}