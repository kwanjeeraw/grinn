#'@import RCurl jsonlite igraph WGCNA RJSONIO
#'@importFrom Hmisc capitalize
#'@importFrom plyr rbind.fill
#'@importFrom stringr str_trim
#'@importFrom reshape2 melt

##list of Cypher to query the path of different relationship types
##new relationship types can be added here
pathList <- list(
  MetaboliteProtein = "UNWIND keyword AS x WITH x MATCH ptw = (from:Metabolite)<-[rel:ASSIGN_TO*]-(to:Protein)",
  ProteinMetabolite = "UNWIND keyword AS x WITH x MATCH ptw = (from:Protein)-[rel:ASSIGN_TO*]->(to:Metabolite)",
  MetaboliteGene = "UNWIND keyword AS x WITH x MATCH ptw = (from:Metabolite)<-[*]-(to:Gene)",
  GeneMetabolite = "UNWIND keyword AS x WITH x MATCH ptw = (from:Gene)-[*]->(to:Metabolite)",
  ProteinGene = "UNWIND keyword AS x WITH x MATCH ptw = (from:Protein)<-[rel:ENCODE*]-(to:Gene)",
  GeneProtein = "UNWIND keyword AS x WITH x MATCH ptw = (from:Gene)-[rel:ENCODE*]->(to:Protein)",
  MetabolitePathway = "UNWIND keyword AS x WITH x MATCH ptw = (from:Metabolite)<-[rel:MAP_TO]-(to:Pathway)",
  PathwayMetabolite = "UNWIND keyword AS x WITH x MATCH ptw = (from:Pathway)-[rel:MAP_TO]->(to:Metabolite)",
  ProteinPathway = "UNWIND keyword AS x WITH x MATCH ptw = (from:Protein)<-[rel:MAP_TO]-(to:Pathway)",
  PathwayProtein = "UNWIND keyword AS x WITH x MATCH ptw = (from:Pathway)-[rel:MAP_TO]->(to:Protein)",
  GenePathway = "UNWIND keyword AS x WITH x MATCH ptw = (from:Gene)<-[rel:MAP_TO]-(to:Pathway)",
  PathwayGene = "UNWIND keyword AS x WITH x MATCH ptw = (from:Pathway)-[rel:MAP_TO]->(to:Gene)"
)

##list of Cypher to query nodes
nodeList <- list(
  exactMatch = "UNWIND keyword AS x WITH x MATCH (node:label) WHERE node.property = x",
  exactCollection = "UNWIND keyword AS x WITH x MATCH (node:label) WHERE ANY(y IN node.property WHERE lower(y) = lower(x))",
  regexMatch = "UNWIND keyword AS x WITH x MATCH (node:label) WHERE lower(str(node.property)) =~ lower(\'.*\'+x+\'.*\')",
  regexCollection = "UNWIND keyword AS x WITH x MATCH (node:label) WHERE ANY(y IN node.property WHERE lower(y) =~ lower(\'.*\'+x+\'.*\'))"
)

##list of property types
##new property can be added here
propertyList <- list(
  stringVal = c("grinn","description","formula","name","inchi","organism","describeIn"),
  listVal = c("xref","synonym","tissue_locations","cellular_locations","biofluid_locations"),
  numericVal = c("averageMass","geneStart","geneEnd")
)

getNodeInfo = function(z, x, y){
  querystring <- paste0("MATCH (node:",z[y],") WHERE node.GID = '",z[x],"' RETURN DISTINCT node")
  result <- curlRequestCypher(querystring)
  if(length(result)>0){
    data.frame(id=result[[1]]$metadata$id,gid=result[[1]]$data$GID,nodename=result[[1]]$data$name,nodetype=result[[1]]$metadata$labels[[1]],xref=paste0(unlist(result[[1]]$data$xref),collapse = "||"), stringsAsFactors = FALSE)
  }else{
    data.frame()
  }
}
getModuleInfo = function(z, x, y){
  querystring <- paste0("MATCH (node:",z[y],") WHERE node.GID = '",z[x],"' RETURN DISTINCT node")
  result <- curlRequestCypher(querystring)
  if(length(result)>0){
    data.frame(id=result[[1]]$metadata$id,gid=result[[1]]$data$GID,nodename=result[[1]]$data$name,nodetype=result[[1]]$metadata$labels[[1]],modulecolor=z["modulecolor"],xref=paste0(unlist(result[[1]]$data$xref),collapse = "||"), stringsAsFactors = FALSE)
  }else{
    data.frame()
  }
}
formatId = function(x, y) {
  ind = which(y$gid == x)
  x = ifelse(length(ind)>0,y$id[ind],x)
}

###for enrichment analysis
#x = combGeneNetwork$edges[1:30,]
#g = graph.edgelist(as.matrix(data.frame(from=unlist(x[,1]),to=unlist(x[,2]))),directed = FALSE)
#g2 = cocitation(g)
#g2[c(2,7,25),c(3,5,8,10)]