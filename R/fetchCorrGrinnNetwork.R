#'Compute a weighted correlation network and expand the network with information from Grinn internal database
#'@description from input omics data e.g. normalized expression data or metabolomics data, it is a one step function to:
#'
#'1. Compute a weighted correlation network of input omics data using WGCNA functions \code{cor} and \code{corPvalueStudent}. 
#'The correlation coefficients are continuous values between -1 (negative correlation) and 1 (positive correlation), with numbers close to 1 or -1, meaning very closely correlated.
#'
#'2. Expand the correlation network using information from the Grinn internal database.
#'The nodes of the correlation network are the keywords input to query the Grinn internal database.
#'The Grinn internal database contains the networks of the following types that can get expanded to: 
#'metabolite-protein, metabolite-protein-gene, metabolite-pathway, protein-gene, protein-pathway and gene-pathway, see also \code{\link{fetchGrinnNetwork}}.
#'@usage fetchCorrGrinnNetwork(datX, datY, corrCoef, pval, method, returnAs, xTo, yTo, filterSource)
#'@param datX data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities. 
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples e.g. normals, tumors. 
#'Require 'nodetype' at the first row to indicate the type of entities in each column. See below for details.
#'@param datY data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities.
#'Use the same format as \code{datX} or it can be NULL. See below for details.
#'@param corrCoef numerical value to define the minimum value of absolute correlation, from 0 to 1, to include edges in the output.
#'@param pval numerical value to define the maximum value of pvalues, to include edges in the output.
#'@param method string to define which correlation is to be used. It can be one of "pearson","kendall","spearman" (default), see \code{\link{cor}}.  
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape", default is "tab". "cytoscape" is the format used in Cytoscape.js
#'@param xTo string of node type. It can be one of "metabolite","protein","gene","pathway". See below for details.
#'@param yTo string of node type. It can be one of "metabolite","protein","gene","pathway". See below for details. 
#'@param filterSource string or list of pathway databases. The argument is required, if \code{xTo} or \code{yTo = "pathway"}.
#'The argument value can be any of "SMPDB","KEGG","REACTOME" or combination of them e.g. list("KEGG","REACTOME").
#'@details datX and datY are matrices in which rows are samples and columns are entities. 
#'
#'If datY is given, then the correlations between the columns of datX and the columns of datY are computed.
#'In this case: 
#'
#'- The correlation network can be expand from datX (by providing a value to \code{xTo}) or datY (by providing a value to \code{yTo}) or both entities to the specified nodetype.
#'
#'Otherwise if datY is not given, the correlations of the columns of datX are computed. 
#'In this case:
#'
#'- The correlation network can be expand from datX entites to a specific nodetype, by providing a value to \code{xTo}.
#'
#'If \code{xTo} or \code{yTo} or both is given, the columns of both datX and datY are required to use grinn ids for extended queries on the Grinn internal database, see \code{\link{convertToGrinnID}} for id conversion.
#'
#'If \code{xTo} = NULL and \code{yTo} = NULL, only the correlation network will be returned.
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references Langfelder P. and Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9:559 
#'@references Dudoit S., Yang YH., Callow MJ. and Speed TP. (2002) Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA, 12:111
#'@references Langfelder P. and Horvath S. Tutorials for the WGCNA package \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html}
#'@export
#'@seealso \code{\link{cor}}, \code{\link{corPvalueStudent}}, \code{\link{fetchGrinnNetwork}}, \url{http://js.cytoscape.org/}
#'@examples
#'# Compute a correlation network of metabolites and expand to a grinn network of metabolite-protein
#'dummy <- rbind(nodetype=rep("metabolite"),t(mtcars))
#'colnames(dummy) <- c('G1.1','G27967','G371','G4.1',paste0('G',sample(400:22000, 28)))
#'result <- fetchCorrGrinnNetwork(datX=dummy, corrCoef=0.7, pval=1e-12, method="spearman", returnAs="tab", xTo="protein")
#'library(igraph)
#'plot(graph.data.frame(result$edges[,1:2], directed=FALSE))
#'# Compute a correlation network of metabolites and proteins and expand to the grinn network of metabolite-pathway and protein-gene
#'dummyX <- rbind(nodetype=rep("metabolite"),t(mtcars)[,1:16])
#'colnames(dummyX) <- c('G1.1','G27967','G371','G4.1',paste0('G',sample(400:22000, 12)))
#'dummyY <- rbind(nodetype=rep("protein"),t(mtcars)[,17:32])
#'colnames(dummyY) <- c('P28845','P08235','Q08AG9','P80365',paste0('P',sample(10000:80000, 12)))
#'result <- fetchCorrGrinnNetwork(datX=dummyX, datY=dummyY, corrCoef=0.7, pval=1e-4, method="spearman", returnAs="json", xTo="pathway", yTo="gene")

fetchCorrGrinnNetwork <- function(datX, datY=NULL, corrCoef=0.5, pval=1e-9, method="spearman", returnAs="tab", 
                                  xTo=NULL, yTo=NULL, filterSource=list()){
  corrnw = fetchCorrNetwork(datX=datX,datY=datY,corrCoef=corrCoef,pval=pval,method=method,returnAs="tab")
  if(nrow(corrnw$nodes)>0){
    nodetypes = tolower(unique(corrnw$nodes$nodetype))
    if(length(nodetypes)>1){#if there are two data types
      if(!is.null(xTo)){
        basicnw4 = fetchGrinnNetwork(txtInput=corrnw$nodes[which(tolower(corrnw$nodes$nodetype)==nodetypes[1]), 1],from=nodetypes[1],
                                   to=xTo,filterSource=filterSource,dbXref="grinn") #relations sources to others
      }else{
        basicnw4 = list(nodes=data.frame(),edges=data.frame())
      }
      if(!is.null(yTo)){
        basicnw5 = fetchGrinnNetwork(txtInput=corrnw$nodes[which(tolower(corrnw$nodes$nodetype)==nodetypes[2]), 1],from=nodetypes[2],
                                     to=yTo,filterSource=filterSource,dbXref="grinn") #relations from targets to others
      }else{
        basicnw5 = list(nodes=data.frame(),edges=data.frame())
      }
      basicnwNodes = rbind(basicnw4$nodes,basicnw5$nodes) #collect all nodes
      basicnwEdges = rbind(basicnw4$edges,basicnw5$edges) #collect all edges
      if(nrow(basicnwEdges)>0){
        basicnwEdges = basicnwEdges[!duplicated(basicnwEdges[,1:ncol(basicnwEdges)]),] #remove duplicated edges
        basicnwNodes = basicnwNodes[!duplicated(basicnwNodes[,1]),] #remove duplicated nodes
        basicnw = list(nodes=basicnwNodes, edges=basicnwEdges)
      }else{
        basicnw = list(nodes=data.frame(),edges=data.frame())
      }
    }else{#if there is only one data type
      if(!is.null(xTo)){
        basicnw = fetchGrinnNetwork(txtInput=corrnw$nodes$id,from=nodetypes,to=xTo,filterSource=filterSource,dbXref="grinn")
      }else{
        basicnw = list(nodes=data.frame(),edges=data.frame())
      }
    }
    #collect node info
    corrattb = data.frame()
    corrattb = plyr::ldply (apply(corrnw$nodes, MARGIN = 1, FUN=getNodeInfo, x = "id", y = "nodetype")) #format nodelist
    corrnw$edges$source = lapply(corrnw$edges$source, FUN=formatId, y = corrattb) #format edgelist
    corrnw$edges$target = lapply(corrnw$edges$target, FUN=formatId, y = corrattb) #format edgelist
    if(nrow(basicnw$nodes)>0){
      cat("Formating and returning combined network ...\n")
      basicnw$edges$corr_coef = 1
      basicnw$edges$pval = 0
      basicnw$edges$direction = 0
      corrnw$edges$relsource = ""
      corrnw$nodes$xref = ""
      corrnw$nodes$gid = corrnw$nodes$id #same ids
      pair = rbind(basicnw$edges,corrnw$edges)
      if(nrow(corrattb)>0){attb = rbind(basicnw$nodes,corrattb,corrnw$nodes)}else{attb = rbind(basicnw$nodes,corrnw$nodes)}
      attb = attb[!duplicated(attb[,2]),]
      cat("Found ",nrow(pair)," relationships...\n")
    }else{#if only correlation network found
      cat("Formating and returning combined network ...\n")
      pair = corrnw$edges
      corrnw$nodes$xref = ""
      corrnw$nodes$gid = corrnw$nodes$id #same ids
      if(nrow(corrattb)>0){attb = rbind(corrattb,corrnw$nodes)}else{attb = corrnw$nodes}
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