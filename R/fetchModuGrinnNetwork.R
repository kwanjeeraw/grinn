#'Compute a network module correlated to a phenotypic feature and expand the network with information from grinn internal database
#'@description from input omics data e.g. normalized expression data or metabolomics data, and phenotypic data it is a one step function to:
#'
#'1. Identify correlation between the input omics data e.g. normalized gene expression data, and phenotypic data e.g. weight.
#'The network module is the output of function \code{fetchWGCNAodule}. The function wraps around important aspects of \pkg{WGCNA} including 
#'\code{blockwiseModules}, \code{cor}, \code{corPvalueStudent}, \code{labeledHeatmap}.
#'These aspects automatically perform correlation network construction, module detection, and display module-phenotype correlations.
#'A module or the combination of modules can be selected from the heatmap of module-phenotype correlations for including in the network output, see more details in \code{\link{fetchWGCNAodule}}.
#'
#'2. Expand the network module using information from Grinn internal database.
#'The nodes of the network module are the keywords input to query the grinn database.
#'Grinn internal database contains the networks of the following types that can get expanded to: 
#'metabolite-protein, metabolite-protein-gene, metabolite-pathway, protein-gene, protein-pathway and gene-pathway, see also \code{\link{fetchGrinnNetwork}}.
#'@usage fetchModuGrinnNetwork(datNorm, datPheno, sfPower, minModuleSize, threshold, returnAs, targetTo, filterSource, organism)
#'@param datNorm data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities. 
#'The column names of datNorm are required to use grinn ids. \code{convertToGrinnID} is provided for id conversion, see \code{\link{convertToGrinnID}}.
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples e.g. normals, tumors. Require 'nodetype' at the first row to indicate the type of entities in each column.
#'@param datPheno data frame containing phenotypic data e.g. weight, age, insulin sensitivity. Columns correspond to phenotypes and rows to samples e.g. normals, tumors. 
#'@param sfPower numerical value of soft-thresholding power for correlation network construction. It is automatically estimated using \code{pickSoftThreshold}, or it can be defined by users.
#'@param minModuleSize numerical value of minimum module size for module detection.
#'@param threshold numerical value to define the minimum value of similarity threshold, from 0 to 1, to include edges in the output.
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape", default is "tab". "cytoscape" is the format used in Cytoscape.js
#'@param targetTo string of node type. It can be one of "metabolite","protein","gene","pathway". By default, the network module will expand to pathways.
#'@param filterSource string or list of pathway databases. The argument is required, if \code{sourceTo} or \code{targetTo = "pathway"}.
#'The argument value can be any of "SMPDB","KEGG","REACTOME" or combination of them e.g. list("KEGG","REACTOME").
#'@param organism string of species in the following format: organism = "'species'". Default is "'Homo sapiens'".
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references Langfelder P. and Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9:559 
#'@references Dudoit S., Yang YH., Callow MJ. and Speed TP. (2002) Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA, 12:111
#'@references Langfelder P. and Horvath S. Tutorials for the WGCNA package \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html}
#'@export
#'@seealso \code{\link{cor}}, \code{\link{corPvalueStudent}}, \code{\link{pickSoftThreshold}}, \code{\link{blockwiseModules}}, \code{\link{labeledHeatmap}}, 
#'\code{\link{exportNetworkToCytoscape}}, \code{\link{fetchWGCNAodule}}, \code{\link{fetchGrinnNetwork}}, \url{http://js.cytoscape.org/}
#'@examples
#'# Compute a correlation of metabolite module to phenotypic data and expand to a grinn network metabolite-pathway, using information from REACTOME and KEGG only
#'library(grinn)
#'data(dummy)
#'data(dummyPheno)
#'result <- fetchModuGrinnNetwork(datNorm=dummy, datPheno=dummyPheno, minModuleSize=5, threshold=0.2, returnAs="tab", targetTo="pathway", filterSource=c("REACTOME","KEGG"))
#'#enter module color(s) seperate by space:yellow brown purple
#'#library(igraph)
#'#plot(graph.data.frame(result$edges[,1:2], directed=FALSE))

fetchModuGrinnNetwork <- function(datNorm, datPheno, sfPower=NULL, minModuleSize = 10, threshold = 0.5, returnAs="tab", 
                                  targetTo, filterSource=list(), organism="'Homo sapiens'"){
  modulenw = fetchWGCNAodule(datNorm=datNorm, datPheno=datPheno, sfPower=sfPower, minModuleSize=minModuleSize, threshold=threshold, returnAs="tab")
  
  if(nrow(modulenw$nodes)>0){
    nodetypes = tolower(unique(modulenw$nodes$nodetype))
    basicnw = fetchGrinnNetwork(txtInput=modulenw$nodes$id,from=nodetypes,to=targetTo,filterSource=filterSource,dbXref="grinn")
    if(nrow(basicnw$nodes)>0){
      cat("Formating and returning combined network ...\n")
      basicnw$edges$corr_coef = 1
      basicnw$nodes$modulecolor = ""
      modulenw$edges$relsource = ""
      modulenw$nodes$xref = ""
      pair = rbind(basicnw$edges,modulenw$edges)
      attb = rbind(basicnw$nodes,modulenw$nodes)
      attb = attb[!duplicated(attb[,1]),]
      colnames(attb) = c("id","nodename","xref","nodetype")
      cat("Found ",nrow(pair)," relationships...\n")
    }else{#if only correlation network found
      cat("Formating and returning combined network ...\n")
      pair = modulenw$edges
      attb = modulenw$nodes
      cat("Found ",nrow(pair)," relationships...\n")
    }
  }else{#if no correlation network found
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