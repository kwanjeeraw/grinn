#' \code{fetchRelation} get path information
#'@description get path information as the output for further uses by \code{formatNetworkOutput}.
#'@seealso \code{formatNetworkOutput}
#'#result <- fetchRelation("http://localhost:7474/db/data/relationship/53")
#'#return start-relation-end

fetchRelation <- function(url){
  out <- tryCatch(
  {
    path = curlRequestUrlToList(url)
    start = curlRequestUrlToList(path$start)
    end = curlRequestUrlToList(path$end)
    type = path$type
    dataSource = path$data$source
    startGID = start$data$GID
    startName = start$data$name
    startXref = paste0(start$data$xref,collapse = "||")
    startLabel = start$metadata$labels[[1]]
    endGID = end$data$GID
    endName = end$data$name
    endXref = paste0(end$data$xref,collapse = "||")
    endLabel = end$metadata$labels[[1]]
    ## Set the name for the class
    relation = list(startGID=startGID, startName=startName, startXref=startXref, startLabel=startLabel, 
                    endGID=endGID, endName=endName, endXref=endXref, endLabel=endLabel, type=type, dataSource=dataSource)
  },
  error=function(e) {
    message(e)
    cat("\n..RETURN empty list of relations")
    out = list() # Choose a return value in case of error
  })    
  return(out)
}

fetchRelation.TRANSACTION <- function(graph){
  out <- tryCatch(
    {
      ## Set the name for the class
      data.frame(t(sapply(graph$graph$relationships, 
                                 function(x) list(source=x$startNode, target=x$endNode, relname=x$type, relsource=x$properties["source"]))))
      
      #relationInfo = list(source=graph$graph$relationships[[1]]$startNode, target=graph$graph$relationships[[1]]$endNode, relname=graph$graph$relationships[[1]]$type, relsource=graph$graph$relationships[[1]]$properties["source"])
    },
    error=function(e) {
      message(e)
      cat("\n..RETURN empty list of relations")
      out = data.frame() # Choose a return value in case of error
    })    
  return(out)
}

fetchNode.TRANSACTION <- function(node){
  out <- tryCatch(
    {
      data.frame(t(sapply(node, 
                          function(x) list(id=x$id, gid=x$properties$GID, nodename=x$properties$name, xref=paste0(x$properties$xref,collapse = "||"), nodetype=x$labels))))
      #nodeInfo = list(id=node$id, gid=node$properties$GID, nodename=node$properties$name, xref=paste0(node$properties$xref,collapse = "||"), nodetype=node$labels)
    },
    error=function(e) {
      message(e)
      cat("\n..RETURN empty list of relations")
      out = data.frame() # Choose a return value in case of error
    })    
  return(out)
}