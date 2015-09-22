#'Set grinn DB location
#'@description set grinn DB location for the currently working environment.
#'@usage setGrinnDb()
#'@param url grinn DB location
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@export
#'@examples
#'setGrinnDb("http://localhost:7474/db/data/cypher")
setGrinnDb <- function(url){
  assign("nld", url, envir = .GlobalEnv)
}