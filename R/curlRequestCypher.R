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
#'#querystring <- "UNWIND ['G7','G8'] AS x WITH x MATCH (node:Metabolite) WHERE lower(node.GID) = lower(x) RETURN DISTINCT node" 
#'#result <- curlRequestCypher(querystring)

curlRequestCypher <- function(querystring){
  h = RCurl::basicTextGatherer()
  tryCatch({
    url = paste0(nld,"cypher")
    RCurl::curlPerform(url=url,
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

curlRequest.TRANSACTION <- function(cypher){
  h = RCurl::basicTextGatherer()
  tryCatch({
    #cypher = "MATCH ptw = (from:Protein)-[:CONTROL]->(to:Gene) WHERE from.GID = 'P14859' RETURN DISTINCT ptw LIMIT 2"
    url = paste0(nld,"transaction/commit")
    #url = "http://localhost:7474/db/data/transaction/commit"
    body = paste0("{\"statements\":[{\"statement\":\"",cypher,"\",\"resultDataContents\":[\"graph\"]}]}")
    RCurl::curlPerform(url=url,
                       userpwd = neu,
                       httpheader = c(Authorization = paste("Basic",RCurl::base64(neu)), 'Content-Type' = "application/json"),
                       postfields=body,
                       writefunction = h$update,
                       verbose = FALSE
    ) 
    result <- RJSONIO::fromJSON(h$value())$results[[1]]$data
  }, error = function(err) {
    message(err)
    result <- list() #return empty if not found
  })
}