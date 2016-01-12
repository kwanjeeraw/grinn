#'Compute a network module correlated to a phenotypic feature and expand the network with information from Grinn internal database
#'@description from input omics data e.g. normalized expression data or metabolomics data, and phenotypic data it is a one step function to:
#'
#'1. Identify correlation between the input omics data e.g. normalized gene expression data, and phenotypic data e.g. weight.
#'The function wraps around important aspects of \pkg{WGCNA} including \code{blockwiseModules}, \code{cor}, \code{corPvalueStudent}, \code{labeledHeatmap}.
#'These aspects automatically perform correlation network construction, module detection, and display module-phenotype correlations.
#'A module or the combination of modules can be selected from the heatmap of module-phenotype correlations for including in the network output, see more details below.
#'
#'2. Expand the network module using information from the Grinn internal database.
#'The nodes of the network module are the keywords input to query the grinn database.
#'The Grinn internal database contains the networks of the following types that can get expanded to: 
#'metabolite-protein, metabolite-protein-gene, metabolite-pathway, protein-gene, protein-pathway and gene-pathway, see also \code{\link{fetchGrinnNetwork}}.
#'@usage fetchModuGrinnNetwork(datX, datPheno, sfPower, minModuleSize, threshold, returnAs, xTo, filterSource)
#'@param datX data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities. 
#'The column names of datX are required to use grinn ids. \code{convertToGrinnID} is provided for id conversion, see \code{\link{convertToGrinnID}}.
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples e.g. normals, tumors. Require 'nodetype' at the first row to indicate the type of entities in each column.
#'@param datPheno data frame containing phenotypic data e.g. weight, age, insulin sensitivity. Columns correspond to phenotypes and rows to samples e.g. normals, tumors. 
#'@param sfPower numerical value of soft-thresholding power for correlation network construction. It is automatically estimated using \code{pickSoftThreshold}, or it can be defined by users.
#'@param minModuleSize numerical value of minimum module size for module detection.
#'@param threshold numerical value to define the minimum value of similarity threshold, from 0 to 1, to include edges in the output.
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape", default is "tab". "cytoscape" is the format used in Cytoscape.js
#'@param xTo string of node type. It can be one of "metabolite","protein","gene","pathway". See below for details.
#'@param filterSource string or list of pathway databases. The argument is required, if \code{sourceTo} or \code{xTo = "pathway"}.
#'The argument value can be any of "SMPDB","KEGG","REACTOME" or combination of them e.g. list("KEGG","REACTOME").
#'@details datX is a matrix in which rows are samples and columns are entities. 
#'The network can be expand from datX entities to the specified nodetype by providing a value to \code{xTo}.
#'
#'If \code{xTo} = NULL, only the network module will be returned.
#'
#'The function encapsulates several methods from \pkg{WGCNA} so that module-phenoty correlation analysis can be fasten. These methods include:
#'
#'- \code{pickSoftThreshold} estimates soft-thresholding powers from scale free topology to build the correlation network.
#'
#'- \code{blockwiseModules} automatically calculates a correlation network and detects modules.
#'Modules are the areas of the network where nodes are densely connected based on their topological overlap measure, see \pkg{WGCNA} for more details. 
#'Each module is labeled by color. By using the color, a module or the combination of modules can be selected ("enter color to the terminal"), for including in the network output.
#'
#'- Module-phenotype correlations and significances are calculated using WGCNA functions \code{cor} and \code{corPvalueStudent}.
#'The correlation coefficients are continuous values between -1 (negative correlation) and 1 (positive correlation), with numbers close to 1 or -1, meaning very closely correlated.
#'
#'- \code{labeledHeatmap} plots a heatmap of module-phenotype correlations. A row is a module and a column is a phynotype. 
#'Each cell presents the corresponding correlation and the pvalue inside parenthesis. Each cell is colored by correlation, red means positive and blue means negative correlation.
#'
#'- \code{exportNetworkToCytoscape} exports a network for using in Cytoscape (\url{http://cytoscape.org/}).
#'The selected module is exported as the network output in which an edge will be included if it similarity threshold above the cutoff, see \code{threshold}.
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references Langfelder P. and Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9:559 
#'@references Dudoit S., Yang YH., Callow MJ. and Speed TP. (2002) Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA, 12:111
#'@references Langfelder P. and Horvath S. Tutorials for the WGCNA package \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html}
#'@export
#'@seealso \code{\link{cor}}, \code{\link{corPvalueStudent}}, \code{\link{pickSoftThreshold}}, \code{\link{blockwiseModules}}, \code{\link{labeledHeatmap}}, 
#'\code{\link{exportNetworkToCytoscape}}, \code{\link{fetchWGCNAModule}}, \code{\link{fetchGrinnNetwork}}, \url{http://js.cytoscape.org/}
#'@examples
#'# Compute a correlation of metabolite module to phenotypic data and expand to a grinn network metabolite-pathway, using information from REACTOME and KEGG only
#'library(grinn)
#'data(dummy)
#'data(dummyPheno)
#'result <- fetchModuGrinnNetwork(datX=dummy, datPheno=dummyPheno, minModuleSize=5, threshold=0.2, returnAs="tab", xTo="pathway", filterSource=c("REACTOME","KEGG"))
#'#enter module color(s) seperate by space:yellow brown purple
#'#library(igraph)
#'#plot(graph.data.frame(result$edges[,1:2], directed=FALSE))

fetchModuGrinnNetwork <- function(datX, datPheno, sfPower=NULL, minModuleSize = 10, threshold = 0.5, returnAs="tab", xTo=NULL, filterSource=list()){
  modulenw = fetchWGCNAModule(datX=datX, datPheno=datPheno, sfPower=sfPower, minModuleSize=minModuleSize, threshold=threshold, returnAs="tab")
  if(nrow(modulenw$nodes)>0){
    nodetypes = tolower(unique(modulenw$nodes$nodetype))
    if(!is.null(xTo)){
      basicnw = fetchGrinnNetwork(txtInput=modulenw$nodes$id,from=nodetypes,to=xTo,filterSource=filterSource,dbXref="grinn")
    }else{
      basicnw = list(nodes=data.frame(),edges=data.frame())
    }
    #collect node info
    moduattb = data.frame()
    moduattb = plyr::ldply (apply(modulenw$nodes, MARGIN = 1, FUN=getModuleInfo, x = "id", y = "nodetype")) #format nodelist
    modulenw$edges$source = lapply(modulenw$edges$source, FUN=formatId, y = moduattb) #format edgelist
    modulenw$edges$target = lapply(modulenw$edges$target, FUN=formatId, y = moduattb) #format edgelist
    if(nrow(basicnw$nodes)>0){
      cat("Formating and returning combined network ...\n")
      basicnw$edges$corr_coef = 1
      basicnw$nodes$modulecolor = ""
      modulenw$edges$relsource = ""
      modulenw$nodes$xref = ""
      modulenw$nodes$gid = modulenw$nodes$id #same ids
      pair = rbind(basicnw$edges,modulenw$edges)
      if(nrow(moduattb)>0){attb = rbind(basicnw$nodes,moduattb,modulenw$nodes)}else{attb = rbind(basicnw$nodes,modulenw$nodes)}
      attb = attb[!duplicated(attb[,2]),]
      cat("Found ",nrow(pair)," relationships...\n")
    }else{#if only correlation network found
      cat("Formating and returning combined network ...\n")
      pair = modulenw$edges
      modulenw$nodes$xref = ""
      modulenw$nodes$gid = modulenw$nodes$id #same ids
      if(nrow(moduattb)>0){attb = rbind(moduattb,modulenw$nodes)}else{attb = modulenw$nodes}
      attb = attb[!duplicated(attb[,2]),]
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