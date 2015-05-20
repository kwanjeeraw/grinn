#' \code{pairToNode} generate an integrated network
#'@description from the list of metabolites, find and connect these metabolites to three type of nodes 
#'including Protein, Gene and Pathway.
#'@usage pairToNode(txtInput, organism, label=c("Protein","Gene","Pathway"), returnAs="json", searchBy="grinn")
#'@param txtInput string containing metabolite ids. Metabolite ids are from the specified database, see \code{searchBy}. Default is grinn id e.g. C00065.
#'Search using grinn ids, txtInput is in the following format: txtInput = "['id'1, 'id2']".
#'Search using other database ids, txtInput must be in the following format: txtInput = "['databaseName:id'1, 'databaseName:id2']".
#'@param organism string of species in the following format: organism = "'species'"
#'@param label character vector of node types to connect to metabolites, default is all three node types: c("Protein","Gene","Pathway")
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be "all", "json" or "cytoscape", but default is "json". "cytoscape" is the format used in Cytoscape.js
#'@param searchBy string of database name. Specify the database of which the ids used as the txtInput, see \code{txtInput}. 
#'It can be "grinn" or "InChI" or "KEGG" or "PubChem" or "ChEBI". Default is "grinn".
#'@return list of nodes and edges encapsulated in json format. Return empty list if found nothing.
#'@note use in grinnWeb only
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references Modify the output generation for Cytoscape.js from \url{https://github.com/cytoscape/r-cytoscape.js/blob/master/cytoscapeJsSimpleNetwork.R}
#'@references Cytoscape.js, Network library for analysis and visualisation \url{http://js.cytoscape.org/}
#'@export
#'@seealso \url{http://js.cytoscape.org/}
#'@examples
#'# Query metabolites by PubChem ids and connect to Protein, Gene and Pathway
#'txtInput <- "['PubChem:3326','PubChem:3436','PubChem:7656']"
#'organism <- "'Homo sapiens'"
#'result <- pairToNode(txtInput, organism, returnAs="all", searchBy="PubChem")
#'# Query metabolites by grinn ids and and connect to Gene and Pathway
#'txtInput <- "['C00024','C00136','C05269']"
#'organism <- "'Homo sapiens'"
#'result <- pairToNode(txtInput, organism, label=c("Gene","Pathway"), returnAs="all", searchBy="grinn")
pairToNode <- function(txtInput, organism, label=c("Protein","Gene","Pathway"), returnAs="json", searchBy="grinn"){
  pair = data.frame() #list of mapped nodes
  attb = data.frame() #list of node attributes
  
  #for each label
  for(i in 1:length(label)){
    #construct query string
    querystring = relationList$pairwise
    if(searchBy == 'grinn'){
      querystring = paste(querystring,'WHERE lower(target.GID) = lower(x) RETURN DISTINCT target.GID, target.name, target.xref, source.GID, source.name, source.xref')
    }else if(searchBy == 'InChI'){
      querystring = paste(querystring,'WHERE lower(target.InChI) = lower(x) RETURN DISTINCT target.GID, target.name, target.xref, source.GID, source.name, source.xref')
    }else{
      querystring = paste(querystring,'WHERE ANY(y IN target.xref WHERE lower(y) = lower(x)) RETURN DISTINCT target.GID, target.name, target.xref, source.GID, source.name, source.xref')
    }
    querystring = gsub("keyword", txtInput, querystring)
    querystring = gsub("species", organism, querystring)
    querystring = gsub("label", label[i], querystring)
    data <- curlRequestCypher(querystring)
    #print(querystring)    
    if(length(data)>0){
      #collect list of mapped nodes and node attributes
      for(j in 1:length(data)){
        pair = rbind(pair,cbind(unlist(data[[j]][1]),unlist(data[[j]][4]),paste("PAIR_TO",label[i],sep = "_"))) #from-to-reltype
        attb = rbind(attb,cbind(unlist(data[[j]][1]),unlist(data[[j]][2]),(data[[j]][3]),"Metabolite")) #attribute from node
        attb = rbind(attb,cbind(unlist(data[[j]][4]),unlist(data[[j]][5]),(data[[j]][6]),label[i])) #attribute to node
      }
      cat("Returning... ",length(data)," pairs to ",label[i], "nodes\n")
    }
    
  }
  
  if(length(pair)>0){
    cat("**..TOTAL of ",nrow(unique(pair))," unique pairs found..**\n")
    colnames(pair) = c("source","target","reltype")
    colnames(attb) = c("GID","nodename","xref","nodetype")
    pair <- unique(pair)
    attb <- unique(attb)
    
    #create cytoscape network
    lsXref = unlist(lapply(attb[,3],FUN=function(x){paste(unlist(x),collapse = "||")}))    
    g <- igraph::graph.edgelist(as.matrix(pair[,-3]), directed = F)  
    g <- igraph::set.edge.attribute(g,name="reltype",value=as.character(pair[,3]))    
    g <- igraph::set.vertex.attribute(g, name="id", value=as.character(attb[,1]))
    g <- igraph::set.vertex.attribute(g, name="nodename", value=as.character(attb[,2]))
    g <- igraph::set.vertex.attribute(g, name="href", value=paste0("node.html?GID=",as.character(attb[,1])))
    g <- igraph::set.vertex.attribute(g, name="xref", value=lsXref)
    g <- igraph::set.vertex.attribute(g, name="nodetype", value=as.character(attb[,4]))
    gudf = igraph::get.data.frame(g,"both")
    rownames(gudf$vertices) = NULL
    names(gudf) = c("nodes","edges")
    colnames(gudf$edges)[1:2] = c("source","target")
    network = createCyNetwork(gudf$nodes, gudf$edges) #return integrated network
    
    #network in json format
    pair = jsonlite::toJSON(pair)
    attb = jsonlite::toJSON(attb)
  }
  else{# if no mapped node found
    print("Returning no data...")
    network = list(nodes="", edges="")
  }
  yield <- switch(returnAs,
                  all = list(nodelist=attb, edgelist=pair, cynetwork=network),
                  json = list(nodelist=attb, edgelist=pair),
                  cytoscape = list(cynetwork=network),
                  stop("data type not included"))
}