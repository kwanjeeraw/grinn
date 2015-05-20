#' \code{integrateNetwork} generate an integrated network
#'@description from the list of metabolites, find and connect these metabolites through four type of relationships 
#'including biochemical reaction, enzyme catalysis, encoding gene and metabolic pathway.
#'@usage integrateNetwork(txtInput, organism, returnAs="json", searchBy="grinn")
#'@param txtInput string containing metabolite ids. Metabolite ids are from the specified database, see \code{searchBy}. Default is grinn id e.g. C00065.
#'Search using grinn ids, txtInput is in the following format: txtInput = "['id'1, 'id2']".
#'Search using other database ids, txtInput must be in the following format: txtInput = "['databaseName:id'1, 'databaseName:id2']".
#'@param organism string of species in the following format: organism = "'species'"
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
#'# Query metabolites by PubChem ids and build a network
#'txtInput <- "['PubChem:3326','PubChem:3436','PubChem:7656']"
#'organism <- "'Homo sapiens'"
#'result <- integrateNetwork(txtInput, organism, returnAs="all", searchBy="PubChem")
#'# Query metabolites by grinn ids and build a network
#'txtInput <- "['C00024','C00136','C05269']"
#'organism <- "'Homo sapiens'"
#'result <- integrateNetwork(txtInput, organism, returnAs="all", searchBy="grinn")
integrateNetwork <- function(txtInput, organism, returnAs="json", searchBy="grinn"){   
  #generate 4 types of networks, add more types if db structure changed
  bcnw = createBiochemNetwork(txtInput,organism,searchBy)
  enznw = connectNodes(txtInput,organism,"enzcatalyze",searchBy)
  genenw = connectNodes(txtInput,organism,"encgene",searchBy)
  ptwnw = connectNodes(txtInput,organism,"pathway",searchBy)
  
  gu = igraph::graph.union(bcnw,enznw,genenw,ptwnw) #merge networks
  
  yield <- tryCatch({
    gudf = igraph::get.data.frame(gu,"both") #get data frame of generated network
    #format node name
    rownames(gudf$vertices) = NULL
    names(gudf) = c("nodes","edges")
    for(i in 1:nrow(gudf$nodes)){
      ind = which(!is.na(gudf$nodes[i,]))[1]
      gudf$nodes$sysname[i] = gudf$nodes[i,ind]
    }
    gudf$nodes = data.frame(id=gudf$nodes$name, name=gudf$nodes$sysname, href=paste0("node.html?GID=",gudf$nodes$name))
    
    #format edge attributes
    gudf$edges$reltype = paste(sep=",",gudf$edges$reltype_1,gudf$edges$reltype_2,
                               gudf$edges$reltype_3,gudf$edges$reltype_4) #more reltype if db structure changed
    gudf$edges = gudf$edges[,c("from","to","biochem","enzcatalyze","encgene","pathway","reltype")] #more intermediates if db structure changed
    colnames(gudf$edges)[1:2] = c("source","target")
    
    #network in cytoscapeJS format
    cynw = createCyNetwork(gudf$nodes, gudf$edges) 
    #network in json format
    pair = jsonlite::toJSON(gudf$edges)
    attb = gudf$nodes
    colnames(attb) = c("GID","nodename","href")
    attb = jsonlite::toJSON(attb)
    switch(returnAs,
           all = list(nodelist=attb, edgelist=pair, cynetwork=cynw),
           json = list(nodelist=attb, edgelist=pair),
           cytoscape = cynw,
           stop("data type not included"))
  }, error = function(err) {
    #on error return empty network
    list(nodes="", edges="")   
  }) # END tryCatch
  
}