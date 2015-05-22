#'Combine a grinn network queried from grinn internal database to a network module correlated to a phenotypic feature
#'@description  from the list of keywords, input omics data e.g. normalized expression data or metabolomics data, and phenotypic data it is a one step function to:
#'
#'1. Build an integrated network (grinn network) by connecting these keywords to a specified node type, see \code{\link{fetchGrinnNetwork}}.
#'The keywords can be any of these node types: metabolite, protein, gene and pathway.
#'Grinn internal database contains the networks of the following types that can be quried: 
#'metabolite-protein, metabolite-protein-gene, metabolite-pathway, protein-gene, protein-pathway and gene-pathway. 
#'
#'2. Identify correlation between the input omics data e.g. normalized gene expression data, and phenotypic data e.g. weight.
#'The network module is the output of function \code{fetchWGCNAodule}. The function wraps around important aspects of \pkg{WGCNA} including 
#'\code{blockwiseModules}, \code{cor}, \code{corPvalueStudent}, \code{labeledHeatmap}.
#'These aspects automatically perform correlation network construction, module detection, and display module-phenotype correlations.
#'A module or the combination of modules can be selected from the heatmap of module-phenotype correlations for including in the network output, see more details in \code{\link{fetchWGCNAodule}}. 
#'
#'3. Combine the grinn network to the network module.
#'@usage fetchGrinnModuNetwork(txtInput, from, to, filterSource, returnAs, dbXref, organism, datNorm, datPheno, sfPower, minModuleSize, threshold)
#'@param txtInput list of keywords containing keyword ids e.g. txtInput = list('id1', 'id2'). 
#'The keyword ids are from the specified database, see \code{dbXref}. Default is grinn id e.g. X371.
#'@param from string of start node. It can be one of "metabolite","protein","gene","pathway".
#'@param to string of end node. It can be one of "metabolite","protein","gene","pathway".
#'@param filterSource string or list of pathway databases. The argument is required, if \code{from} or \code{to = "pathway"}, see \code{from} and \code{to}.
#'The argument value can be any of "SMPDB","KEGG","REACTOME" or combination of them e.g. list("KEGG","REACTOME").  
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape", default is "tab". "cytoscape" is the format used in Cytoscape.js
#'@param dbXref string of database name. Specify the database name used for the txtInput ids, see \code{txtInput}. 
#'It can be one of "grinn","chebi","kegg","pubchem","inchi","hmdb","smpdb","reactome","uniprot","ensembl","entrezgene". Default is "grinn".
#'If pubchem is used, it has to be pubchem SID (substance ID).
#'@param organism string of species in the following format: organism = "'species'". Default is "'Homo sapiens'".
#'@param datNorm data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities. 
#'The column names of datNorm are required to use grinn ids. \code{convertToGrinnID} is provided for id conversion, see \code{\link{convertToGrinnID}}.
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples e.g. normals, tumors. Require 'nodetype' at the first row to indicate the type of entities in each column.
#'@param datPheno data frame containing phenotypic data e.g. weight, age, insulin sensitivity. Columns correspond to phenotypes and rows to samples e.g. normals, tumors. 
#'@param sfPower numerical value of soft-thresholding power for correlation network construction. It is automatically estimated using \code{pickSoftThreshold}, or it can be defined by users.
#'@param minModuleSize numerical value of minimum module size for module detection.
#'@param threshold numerical value to define the minimum value of similarity threshold, from 0 to 1, to include edges in the output.
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references Langfelder P. and Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9:559 
#'@references Dudoit S., Yang YH., Callow MJ. and Speed TP. (2002) Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA, 12:111
#'@references Langfelder P. and Horvath S. Tutorials for the WGCNA package \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html}
#'@export
#'@seealso \code{\link{cor}}, \code{\link{corPvalueStudent}}, \code{\link{pickSoftThreshold}}, \code{\link{blockwiseModules}}, \code{\link{labeledHeatmap}}, 
#'\code{\link{exportNetworkToCytoscape}}, \code{\link{fetchWGCNAodule}}, \code{\link{fetchGrinnNetwork}}, \url{http://js.cytoscape.org/}
#'@examples
#'# Create metabolite-gene network from the list of metabolites using grinn ids and combine the grinn network to a correlation of metabolite module to phenotypic data
#'kw <- c('X160','X300','X371','X16414','X17191')
#'library(grinn)
#'data(dummy)
#'data(dummyPheno)
#'result <- fetchGrinnModuNetwork(txtInput=kw, from="metabolite", to="gene", datNorm=dummy, datPheno=dummyPheno, minModuleSize=5, threshold=0.2)
#'# enter module color(s) seperate by space:yellow brown purple
#'library(igraph)
#'plot(graph.data.frame(result$edges[,1:2], directed=FALSE))

fetchGrinnModuNetwork <- function(txtInput, from, to, filterSource=list(), returnAs="tab", dbXref="grinn", organism="'Homo sapiens'", 
                                  datNorm, datPheno, sfPower=NULL, minModuleSize = 10, threshold = 0.5){
  basicnw = fetchGrinnNetwork(txtInput=txtInput,from=from,to=to,filterSource=filterSource,returnAs=returnAs,dbXref=dbXref,organism=organism)
  modulenw = fetchWGCNAodule(datNorm=datNorm, datPheno=datPheno, sfPower=sfPower, minModuleSize=minModuleSize, threshold=threshold, returnAs=returnAs)
  if(nrow(basicnw$nodes)>0 && nrow(modulenw$nodes)>0){
    cat("Formating and returning combined network ...\n")
    basicnw$edges$corr_coef = 1
    basicnw$nodes$modulecolor = ""
    modulenw$edges$relsource = ""
    modulenw$nodes$xref = ""
    pair = rbind(basicnw$edges,modulenw$edges)
    attb = rbind(basicnw$nodes,modulenw$nodes)
    attb = attb[!duplicated(attb[,1]),]
    colnames(attb) = c("id","nodename","xref","nodetype","modulecolor")
    cat("Found ",nrow(pair)," relationships...\n")
  }else if(nrow(basicnw$nodes)>0 && nrow(modulenw$nodes)==0){
    cat("Formating and returning combined network ...\n")
    pair = basicnw$edges
    attb = basicnw$nodes
    cat("Found ",nrow(pair)," relationships...\n")
  }else if(nrow(basicnw$nodes)==0 && nrow(modulenw$nodes)>0){
    cat("Formating and returning combined network ...\n")
    pair = modulenw$edges
    attb = modulenw$nodes
    cat("Found ",nrow(pair)," relationships...\n")
  }else{# if no mapped node found
    print("Returning no data...")
    pair = data.frame()
    attb = data.frame()
    cynetwork = list(nodes="", edges="")
  }
  out = switch(returnAs,
               tab = list(nodes=attb, edges=pair),
               json = list(nodes=jsonlite::toJSON(attb), edges=jsonlite::toJSON(pair)),
               cytoscape = createCyNetwork(attb, pair),
               stop("incorrect return type"))
}