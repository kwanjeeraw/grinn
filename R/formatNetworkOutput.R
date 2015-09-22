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
  cat("Formating and returning queried network ...\n")
  if(length(network)<6000){
    for(i in 1:length(network)){   
      path = lapply(unlist(network[[i]]$relationships),fetchRelation)
      tmppair = data.frame(t(sapply(path,c)))
      reltype = paste0(tmppair[,4],"_",tmppair[,8])
      tmppair = unname(tmppair)
      attb = rbind(attb,as.matrix(tmppair[,c(1:4)]),as.matrix(tmppair[,c(5:8)]))
      tmppair = tmppair[,c(1,5,10,9)]
      pair = rbind(pair,cbind(tmppair,reltype)) #from-to-relname-reltype 
    }
  }else{
    cat("Found ",length(network)," but returning 6000 of all relationships...\n")
    for(i in 1:6000){   
      path = lapply(unlist(network[[i]]$relationships),fetchRelation)
      tmppair = data.frame(t(sapply(path,c)))
      reltype = paste0(tmppair[,4],"_",tmppair[,8])
      tmppair = unname(tmppair)
      attb = rbind(attb,as.matrix(tmppair[,c(1:4)]),as.matrix(tmppair[,c(5:8)]))
      tmppair = tmppair[,c(1,5,10,9)]
      pair = rbind(pair,cbind(tmppair,reltype)) #from-to-relname-reltype 
    }
  }
  if(length(pair)>0){
    pair <- unique(pair)
    attb <- unique(attb)
    cat("Found ",nrow(pair)," relationships...\n")
    colnames(pair) = c("source","target","relsource","relname","reltype")
    colnames(attb) = c("id","nodename","xref","nodetype") 
    #network in cytoscapeJS format
    #cynetwork = createCyNetwork(attb, pair) 

  }else{# if no mapped node found
    print("Returning no data...")
    cynetwork = list(nodes="", edges="")
  }

  out = switch(returnAs,
                  tab = list(nodes=attb, edges=pair),
                  json = list(nodes=jsonlite::toJSON(attb), edges=jsonlite::toJSON(pair)),
                  cytoscape = createCyNetwork(attb, pair),
                  stop("incorrect return type"))
}