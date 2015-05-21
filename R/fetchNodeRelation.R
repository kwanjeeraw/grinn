#' \code{fetchNodeRelation} format get node relationships
#'@description get node relationships as the output for further uses by \code{formatNodeOutput}.
#'@seealso \code{formatNodeOutput}
#'#result <- fetchNodeRelation("http://localhost:7474/db/data/node/53/relationships/in")
#'#return start-relation-end

fetchNodeRelation <- function(url){
  path = curlRequestUrlToList(url)
  relations = list()
  if(length(path)>0){
    for(i in 1:length(path)){
      start = curlRequestUrlToList(path[[i]]$start)
      end = curlRequestUrlToList(path[[i]]$end)
      type = path[[i]]$type
      dataSource = path[[i]]$data$source
      startGID = start$data$GID
      startName = start$data$name
      startXref = paste0(start$data$xref,collapse = "||")
      startLabel = start$metadata$labels[[1]]
      endGID = end$data$GID
      endName = end$data$name
      endXref = paste0(end$data$xref,collapse = "||")
      endLabel = end$metadata$labels[[1]]
      ## Set the name for the class
      relation = data.frame(startGID=startGID, startName=startName, startXref=startXref, startLabel=startLabel, 
                      endGID=endGID, endName=endName, endXref=endXref, endLabel=endLabel, type=type, dataSource=dataSource)
      relations = rbind(relations,relation)
    }
    relations <- unique(relations)
  }else{
    #cat("\n..RETURN empty list")
    relations = list() # Choose a return value in case of error
  } 
  return(relations)
}