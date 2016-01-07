##@import RCurl jsonlite igraph opencpu
##use in grinnWeb only
##list of Cypher to query the path of different relationship types
relationList <- list(
  biochem = "MATCH (ptw:Pathway{organism:species})-[:HAS]->(rx:Reaction) WITH rx MATCH left<-[:TRANSFORM]-(rx)-[:PRODUCE]->right",
  enzcatalyze = "UNWIND keyword AS x WITH x MATCH ptw = allShortestPaths((target:Metabolite)<-[*]-(source:Protein {organism:species}))",
  encgene = "UNWIND keyword AS x WITH x MATCH ptw = allShortestPaths((target:Metabolite)<-[*]-(source:Gene {organism:species}))",
  pathway = "UNWIND keyword AS x WITH x MATCH ptw = allShortestPaths((target:Metabolite)<-[*]-(source:Pathway {organism:species}))",
  pairwise = "UNWIND keyword AS x WITH x MATCH ptw = allShortestPaths((target:Metabolite)<-[*]-(source:label {organism:species}))"
)

##function to get complete node information
getNodeInfo <- function(data){
  for(i in 1:length(data)){
    outurl = data[[i]]$outgoing_relationships
    data[[i]] = RCurl::merge.list(data[[i]], lapply(getRelInfo(outurl,"out"),FUN=list))
    inurl = data[[i]]$incoming_relationships
    data[[i]] = RCurl::merge.list(data[[i]], lapply(getRelInfo(inurl,"in"),FUN=list))
  }
  result <- data #output node info + its relationships info
}

##function to get all relationship information of a node
getRelInfo <- function(url,direction){
  result <- curlRequestUrlToList(url)
  if(direction=="out"){#outgoing rel
    rellist.names = c("end_name","end_GID","out_type","out_data") #set names
    rel <- sapply(rellist.names,function(x) NULL) #set default relationship value to null
    if(length(result)>0){
      for(i in 1:length(result)){
        dirurl = result[[i]]$end
        rel$end_name = c(rel$end_name, curlRequestUrlToList(dirurl)$data$name)
        rel$end_GID = c(rel$end_GID, curlRequestUrlToList(dirurl)$data$GID)
        rel$out_type = c(rel$out_type, result[[i]]$type)
        
        if(length(result[[i]]$data)>0){
          rel$out_data = c(rel$out_data, jsonlite::toJSON(result[[i]]$data))
        }
        else{
          rel$out_data = c(rel$out_data, jsonlite::toJSON(""))#no property default value
        }
      }
    }
  }#out
  else{#incoming rel
    rellist.names = c("start_name","start_GID","in_type","in_data") #set names
    rel <- sapply(rellist.names,function(x) NULL) #set default relationship value to null
    if(length(result)>0){
      for(i in 1:length(result)){
        dirurl = result[[i]]$start
        rel$start_name = c(rel$start_name, curlRequestUrlToList(dirurl)$data$name)
        rel$start_GID = c(rel$start_GID, curlRequestUrlToList(dirurl)$data$GID)
        rel$in_type = c(rel$in_type, result[[i]]$type)
        
        if(length(result[[i]]$data)>0){
          rel$in_data = c(rel$in_data, jsonlite::toJSON(result[[i]]$data))
        }
        else{
          rel$in_data = c(rel$in_data, jsonlite::toJSON(""))#no property default value
        }
      }
    }
  }#in
  rel #output relationships info
}

##function to create metabolite - RX - metabolite network
createBiochemNetwork <- function(txtInput, organism, searchBy){ 
  #construct query string
  querystring = relationList$biochem
  if(searchBy == 'grinn'){
    querystring = paste(querystring,'WHERE ANY(y IN keyword WHERE lower(y) = lower(left.GID)) AND ANY(y IN keyword WHERE lower(y) = lower(right.GID)) RETURN left.GID, left.name, right.GID, right.name, rx.GID, rx.name ORDER BY left.GID')
  }else if(searchBy == 'InChI'){
    querystring = paste(querystring,'WHERE ANY(y IN keyword WHERE lower(y) = lower(left.InChI)) AND ANY(y IN keyword WHERE lower(y) = lower(right.InChI)) RETURN left.GID, left.name, right.GID, right.name, rx.GID, rx.name ORDER BY left.GID')
  }else{
    querystring = paste(querystring,'WHERE ANY(y IN keyword WHERE ANY(z IN left.xref WHERE lower(z) = lower(y))) AND ANY(y IN keyword WHERE ANY(z IN right.xref WHERE lower(z) = lower(y))) RETURN left.GID, left.name, right.GID, right.name, rx.GID, rx.name ORDER BY left.GID')
  }
  querystring = gsub("keyword", txtInput, querystring)
  querystring = gsub("species", organism, querystring)
  
  # querystring = "MATCH (ptw:Pathway{organism:\"Homo sapiens\"})-[:HAS]->(rx:Reaction) WITH rx MATCH left<-[:TRANSFORM]-(rx)-[:PRODUCE]->right WHERE ANY(y IN [\"C00024\",\"C00136\",\"C00010\",\"C05269\"] WHERE lower(y) = lower(left.GID)) AND ANY(y IN [\"C00024\",\"C00136\",\"C00010\",\"C05269\"] WHERE lower(y) = lower(right.GID)) RETURN left.GID, left.name, right.GID, right.name, rx.GID, rx.name ORDER BY left.GID"
  #print(querystring)
  data <- curlRequestCypher(querystring) #table of left.GID, left.name, right.GID, right.name, rx.GID, rx.name
  
  result <- tryCatch({   
    #create a temp graph from query result
    g <- igraph::graph.edgelist(matrix(data[,c(1,3)], ncol=2), directed = F)
    #set node name
    nodes = unique(rbind(data[,1:2], data[,3:4]))
    findNodeName = function(x){ which(nodes == x) } #internal function to find node name
    ind = lapply(igraph::V(g)$name,findNodeName)
    igraph::V(g)$sysname = nodes[unlist(ind),2]
    #set edge attribute
    igraph::E(g)$biochem = paste0(data[,5],'|',data[,6])
    #create network
    g2 = igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list)
    nw = delete.vertices(g2, v = which(degree(g2,V(g2))==0)) #remove node with no degree
    if(igraph::ecount(nw)>0){# if graph is not empty
      igraph::E(nw)$biochem = lapply(igraph::E(nw)$biochem,unique) #remove duplicate
      combineAttb = function(x){ paste0(unlist(unlist(x)),collapse="||") } #combine reactions, format: GID|name||GID|name, each reaction seperated by ||
      igraph::E(nw)$biochem = lapply(igraph::E(nw)$biochem,combineAttb)
      cat("Returning... ",length(igraph::E(nw))," biochem relationships\n")
    }
    else{# if graph is empty
      igraph::E(nw)$biochem = ""
    }
    igraph::E(nw)$reltype = "biochem" #set relationship type as edge attribute
    result <- nw #return biochem network
  }, error = function(err) {
    #on error return empty network
    result <- igraph::graph.empty(n=0, directed=FALSE)
    igraph::V(result)$name = ""
    igraph::V(result)$sysname = ""
    result <- igraph::set.edge.attribute(result, name= "biochem", value= "")
    result <- igraph::set.edge.attribute(result, name= "reltype", value= "biochem")
  }) # END tryCatch
  #..code to create Cy network..#
  #   df = igraph::get.data.frame(nw,"both") #get data frame of generated network
  #   #format name
  #   rownames(df$vertices) = NULL
  #   names(df) = c("nodes","edges")
  #   names(df$nodes) = c("id","name")
  #   colnames(df$edges)[1:2] = c("source","target")
  #   network <- createCyNetwork(df$nodes, df$edges) 
}

##function to create metabolite - PRT|GN|PTW - metabolite network
connectNodes <- function(txtInput, organism, reltype, searchBy){
  #construct query string
  querystring = relationList[reltype]
  if(searchBy == 'grinn'){
    querystring = paste(querystring,'WHERE lower(target.GID) = lower(x) RETURN DISTINCT target.GID, target.name, source.GID, source.name ORDER by source.GID')
  }else if(searchBy == 'InChI'){
    querystring = paste(querystring,'WHERE lower(target.InChI) = lower(x) RETURN DISTINCT target.GID, target.name, source.GID, source.name ORDER by source.GID')
  }else{
    querystring = paste(querystring,'WHERE ANY(y IN target.xref WHERE lower(y) = lower(x)) RETURN DISTINCT target.GID, target.name, source.GID, source.name ORDER by source.GID')
  }
  querystring = gsub("keyword", txtInput, querystring)
  querystring = gsub("species", organism, querystring)
  
  # querystring = "UNWIND [\"met5\",\"met1\",\"met8\"] AS x WITH x MATCH (source:Protein{organism:\"Homo sapiens\"}), (target:Metabolite) WHERE lower(target.GID) = lower(x) 
  # WITH target, source MATCH ptw = target<-[:TRANSFORM|PRODUCE]-()<-[:CATALYZE]-source 
  #   RETURN target.GID, target.name, source.GID, source.name ORDER by source.GID"
  #print(querystring)
  data <- curlRequestCypher(querystring)
  result <- tryCatch({ 
    if(typeof(data)=="list"){data = t(data.frame(data))} #check type of data, change to dataframe of characters if need
    data <- data[!duplicated(data), ] #table of target.GID, target.name, source.GID, source.name
    #create a temp graph from query result
    g <- igraph::graph.edgelist(matrix(data[,c(1,3)], ncol=2), directed = F)
    toNodes = unique(data[,1])
    
    #find if a node connects to a node from cocitation, connect nodes and create a network
    cocite = igraph::cocitation(g)[toNodes,toNodes] #get adj. matrix of nodes from cocitation
    nw <- igraph::graph.adjacency(cocite, mode = "undirected")
    nw <- igraph::simplify(nw, remove.multiple = T, remove.loops = T)
    el = igraph::get.edgelist(nw) #remove unconnected nodes
    nw <- igraph::graph.edgelist(el, directed = F)
    #set node name
    findNodeName = function(x){ which(data == x)[1] } #internal function to find node name
    ind = lapply(igraph::V(nw)$name,findNodeName)
    igraph::V(nw)$sysname = data[unlist(ind),2]
    #get intermediates
    gr = by(data,data[,3],data.frame)
    findRow = function(x){ which(nrow(x)>1) } #internal function to screen out path with no-intermediate node
    ind = lapply(gr, findRow)
    intmd = gr[names(unlist(ind))]
    findIntermediateNode = function(x){ #internal function to collect intermediates
      lsImd = ""
      for (i in 1:length(intmd)){
        if(all(x %in% intmd[i][[1]]$V1)){
          #format: GID|name||GID|name, each intermediate seperated by ||
          tmp = paste0(as.character(intmd[i][[1]]$V3[1]),'|',as.character(intmd[i][[1]]$V4[1]))
          lsImd = paste0(lsImd,tmp,sep="||")
          #lsImd = c(lsImd, paste0(as.character(intmd[i][[1]]$V3[1]),'|',as.character(intmd[i][[1]]$V4[1])))
        }  
      }
      lsImd = substr(lsImd, 1, nchar(lsImd)-2)
    }
    
    df = igraph::get.data.frame(nw) #get data frame of generated network
    nw = igraph::set.edge.attribute(nw, name= reltype, value= apply(df,1,findIntermediateNode)) #set intermediates as edge attribute, throw error if unfound
    igraph::E(nw)$reltype = reltype #set relationship type as edge attribute
    cat("Returning... ",length(igraph::E(nw))," ",reltype," relationships\n")
    result <- nw #return connected nodes
  }, error = function(err) {
    #on error return empty network
    result <- igraph::graph.empty(n=0, directed=FALSE)
    igraph::V(result)$name = ""
    igraph::V(result)$sysname = ""
    result <- igraph::set.edge.attribute(result, name= reltype, value= "")
    result <- igraph::set.edge.attribute(result, name= "reltype", value= reltype)
  }) # END tryCatch
  #..code to create Cy network..#  
  #   df = igraph::get.data.frame(nw,"both") #get data frame of generated network
  #   #format name
  #   rownames(df$vertices) = NULL
  #   names(df) = c("nodes","edges")
  #   names(df$nodes) = c("id","name")
  #   colnames(df$edges)[1:2] = c("source","target")
  #   network <- createCyNetwork(df$nodes, df$edges)
} 