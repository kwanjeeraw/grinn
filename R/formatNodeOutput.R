#' \code{formatNodeOutput} format resulting node
#'@description Format node output for further uses by \code{fetchNode}. Retrieve also node relationships using \code{fetchNodeRelation}.
#'@seealso \code{fetchNode}, \code{fetchNodeRelation}
#'#querystring="UNWIND ['G371','G783'] AS x WITH x MATCH (node:Metabolite) WHERE lower(node.GID) = lower(x) RETURN DISTINCT node"
#'#node = curlRequestCypher(querystring)
#'#formatNodeOutput(node,returnAs="list")

formatNodeOutput <- function(node,returnAs){
  nodeObj = list(node=list(),incoming_relationships=list(),outgoing_relationships=list())
  nodeObjJson = list(node=list(),incoming_relationships=list(),outgoing_relationships=list())
  nodes = list()
  nodesJson = list()
  for(i in 1:length(node)){
    #list
    nodeObj$node = node[[i]]$data
    nodeObj$incoming_relationships = fetchNodeRelation(node[[i]]$incoming_relationships)
    nodeObj$outgoing_relationships = fetchNodeRelation(node[[i]]$outgoing_relationships)
    nodes = append(nodes,list(nodeObj))
    #json
    nodeObjJson$node = jsonlite::toJSON(nodeObj$node)
    nodeObjJson$incoming_relationships = jsonlite::toJSON(nodeObj$incoming_relationships)
    nodeObjJson$outgoing_relationships = jsonlite::toJSON(nodeObj$outgoing_relationships)
    nodesJson = append(nodesJson,list(nodeObjJson))
  }
  out = switch(returnAs,
               list = nodes,
               json = nodesJson,
               stop("return type not included"))
}