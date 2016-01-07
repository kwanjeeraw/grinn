#' \code{createCyNetwork} format resulting network to cytoscapeJS format
#'@description Modify the R function from r-cytoscape.js, to return json output for generating Cytoscape.js network 
#'@references Modify the output generation for Cytoscape.js from \url{https://github.com/cytoscape/r-cytoscape.js/blob/master/cytoscapeJsSimpleNetwork.R}
#'@seealso \url{https://github.com/cytoscape/r-cytoscape.js/blob/master/cytoscapeJsSimpleNetwork.R}

createCyNetwork <- function(nodeData, edgeData, 
                            nodeColor="#B0B0B0", nodeShape="ellipse") {  
  
  # There must be nodes and nodeData must have at least id and nodename columns
  if(nrow(nodeData) == 0 || !(all(c("id", "nodename") %in% names(nodeData)))) {
    return(list(nodes="", edges=""))
  }
  # There must be edges and edgeData must have at least source and target columns
  if(nrow(edgeData) == 0 || !(all(c("source", "target") %in% names(edgeData)))) {
    return(list(nodes="", edges=""))
  }
  
  # NODES
  ## Add color/shape columns if not present
  if(!("color" %in% colnames(nodeData))) {
    nodeData$color <- rep(nodeColor, nrow(nodeData))
  }
  if(!("shape" %in% colnames(nodeData))) {
    nodeData$shape <- rep(nodeShape, nrow(nodeData))
  }
  nodeEntries <- NULL
  for(i in 1:nrow(nodeData)) {   
    tmpEntries <- NULL
    for(col in colnames(nodeData)) {
      tmp2 <- paste0("\"", col, "\":\"", nodeData[i, col], "\"")
      tmpEntries <- c(tmpEntries, tmp2)
    }
    tmpEntries <- paste(tmpEntries, collapse=", ")
    tmp <- paste0("{ \"data\": { ", tmpEntries, "} }")  
    nodeEntries <- c(nodeEntries, tmp)
  }
  nodeEntries <- paste(nodeEntries, collapse=", ")
  
  # EDGES 
  edgeEntries <- NULL
  for(i in 1:nrow(edgeData)) {   
    tmpEntries <- NULL
    for(col in colnames(edgeData)) {
      tmp2 <- paste0("\"", col, "\":\"", edgeData[i, col], "\"")
      tmpEntries <- c(tmpEntries, tmp2)
    }
    tmpEntries <- paste(tmpEntries, collapse=", ")
    tmp <- paste0("{ \"data\": { ", tmpEntries, "} }")
    edgeEntries <- c(edgeEntries, tmp)
  }
  edgeEntries <- paste(edgeEntries, collapse=", ")
  network <- list(nodes=nodeEntries, edges=edgeEntries)
  
  #print(network)
  return(network)
}