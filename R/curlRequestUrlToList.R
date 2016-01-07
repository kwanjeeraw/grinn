#'Execute query from URI
#'@description execute database query using a URI input. The function calls curlPerform from \pkg{RCurl} for sending HTTP request to the graph database.
#'@usage curlRequestUrlToList(url)
#'@param url string of URI
#'@return list of quried node or relationships. Return \code{NULL} if found nothing.
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references The RCurl package \url{http://www.inside-r.org/packages/cran/RCurl}
#'@seealso \code{\link{curlRequestCypher}}, \code{\link{curlRequestUrlToDF}}, \code{\link{curlPerform}}, \url{http://neo4j.com/docs/milestone/introduction.html}
#'@examples
#'# Get the list of relationships of node 21 
#'#url <- "http://localhost:7474/db/data/node/21/relationships/all" 
#'#result <- curlRequestUrlToList(url)

curlRequestUrlToList <- function(url) {
  h = RCurl::basicTextGatherer()
  tryCatch({
    RCurl::curlPerform(url=url,
                       userpwd = neu,
                       httpheader = c(Authorization = paste("Basic",RCurl::base64(neu))),
                       writefunction = h$update,
                       verbose = FALSE
    ) 
    result <- jsonlite::fromJSON(h$value(), simplifyVector = FALSE) #return as list
  }, error = function(err) {
    print(err)
    result <- NULL #return NULL if not found
  }) # END tryCatch
}