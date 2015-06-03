#'Execute Cypher query
#'@description execute query using Cypher as an input. The function calls curlPerform from \pkg{RCurl} for sending HTTP request to the graph database.
#'@usage curlRequestCypher(querystring)
#'@param querystring string of Cypher
#'@return list of quried results. Return \code{NULL} if found nothing.
#'@note Need to set db location
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references The RCurl package \url{http://www.inside-r.org/packages/cran/RCurl}
#'@seealso \code{\link{curlRequestUrlToDF}}, \code{\link{curlRequestUrlToList}}, \code{\link{curlPerform}}, \url{http://neo4j.com/docs/milestone/introduction.html}
#'@examples
#'# Query metabolites by database id 
#'#querystring <- "UNWIND ['X7','X8'] AS x WITH x MATCH (node:Metabolite) WHERE lower(node.GID) = lower(x) RETURN DISTINCT node" 
#'#result <- curlRequestCypher(querystring)

curlRequestCypher <- function(querystring){
  h = RCurl::basicTextGatherer()
  tryCatch({
    RCurl::curlPerform(url=nld,
                       userpwd = neu,
                       httpheader = c(Authorization = paste("Basic",RCurl::base64(neu))),
                       postfields=paste('query',RCurl::curlEscape(querystring), sep='='),
                       writefunction = h$update,
                       verbose = FALSE
    ) 
    data <- jsonlite::fromJSON(h$value()) 
    result <- data$data #return as data.frame
  }, error = function(err) {
    print(err)
    result <- NULL #return NULL if not found
  }) # END tryCatch
}