#' \code{formatNetworkOutput} format resulting network
#'@description Format network output for further uses by \code{fetchGrinnNetwork}. Retrieve also node relationships using \code{fetchNodeRelation}.
#'@seealso \code{fetchGrinnNetwork}, \code{fetchNodeRelation}
#'#querystring="UNWIND ['G371','G783'] AS x WITH x MATCH ptw = (from:Metabolite)<-[rel:MAP_TO]-(to:Pathway) 
#'#WHERE lower(from.GID) = lower(x) RETURN DISTINCT ptw"
#'#network = curlRequestCypher(querystring)
#'#formatNetworkOutput(network,returnAs="tab")

formatNetworkOutput <- function(network,returnAs){
  pair = data.frame() #list of mapped nodes
  attb = data.frame() #list of node attributes
  if(length(network)<30000){
    cat("Formating and returning network of approximate size", length(network),"...\n")
    pair = plyr::rbind.fill(lapply(network, function(x) fetchRelation.TRANSACTION(x)))
    attb = plyr::rbind.fill(lapply(network, function(x) fetchNode.TRANSACTION(x$graph$nodes)))
    pair = unique(pair)
    attb = unique(attb)
  }else{
    cat("Found ",length(network)," but returning network of size 30000...\n")
    pair = plyr::rbind.fill(lapply(network[1:30000], function(x) fetchRelation.TRANSACTION(x)))
    attb = plyr::rbind.fill(lapply(network[1:30000], function(x) fetchNode.TRANSACTION(x$graph$nodes)))
    pair = unique(pair)
    attb = unique(attb)
  }

  out = switch(returnAs,
                  tab = list(nodes=attb, edges=pair),
                  json = list(nodes=jsonlite::toJSON(attb), edges=jsonlite::toJSON(pair)),
                  cytoscape = createCyNetwork(attb, pair),
                  stop("incorrect return type"))
}