#'Compute a differential correlation network and expand the network with information from Grinn internal database
#'@description from input omics data e.g. normalized expression data or metabolomics data, it is a one step function to:
#'
#'1. Compute a differential correlation network of input omics data from two conditions, see \code{datX1}, \code{datX2}, \code{datY1}, \code{datY2}.
#'Correlation coefficients, pvalues and relation directions among entities in each condition are calculated using WGCNA functions \code{cor} and \code{corPvalueStudent}.
#'The correlation coefficients are continuous values between -1 (negative correlation) and 1 (positive correlation), with numbers close to 1 or -1, meaning very closely correlated.
#'Then correlation coefficients are test for differential correlations using Fisher's z-test based on \pkg{DiffCorr}.
#'
#'2. Expand the differential correlation network using information from the Grinn internal database.
#'The nodes of the differential correlation network are the keywords input to query the grinn database.
#'The Grinn internal database contains the networks of the following types that can get expanded to: 
#'metabolite-protein, metabolite-protein-gene, metabolite-pathway, protein-gene, protein-pathway and gene-pathway, see also \code{\link{fetchGrinnNetwork}}.
#'@usage fetchDiffCorrGrinnNetwork(datX1,  datX2, datY1, datY2, pDiff, method, returnAs, xTo, yTo, filterSource)
#'@param datX1 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of one condition. 
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples. 
#'Require 'nodetype' at the first row to indicate the type of entities in each column. See below for details.
#'@param datX2 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of another condition. 
#'Use the same format as \code{datX1}.
#'@param datY1 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of one condition.
#'Use the same format as \code{datX1}. If there is only one type of dataset, it can be NULL. See below for details.
#'@param datY2 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of another condition.
#'Use the same format as \code{datX1}. If there is only one type of dataset, it can be NULL. See below for details.
#'@param pDiff numerical value to define the maximum value of pvalues (pvalDiff), to include edges in the output.
#'@param method string to define which correlation is to be used. It can be one of "pearson","kendall","spearman", see \code{\link{cor}}.  
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape", default is "tab". "cytoscape" is the format used in Cytoscape.js
#'@param xTo string of node type. It can be one of "metabolite","protein","gene","pathway". See below for details.
#'@param yTo string of node type. It can be one of "metabolite","protein","gene","pathway". See below for details.
#'@param filterSource string or list of pathway databases. The argument is required, if \code{xTo} or \code{yTo = "pathway"}.
#'The argument value can be any of "SMPDB","KEGG","REACTOME" or combination of them e.g. list("KEGG","REACTOME").
#'@details To calculate the differential correlation network, require the input data from two conditions; 1 and 2. 
#'The input data are matrices in which rows are samples and columns are entities.
#'For each condition, if datY is given, then the correlations between the columns of datX and the columns of datY are computed before testing.
#'In this case: 
#'
#'- The differential correlation network can be expand from datX (by providing a value to \code{xTo}) or datY (by providing a value to \code{yTo}) or both entities to the specified nodetype.
#'
#'Otherwise if datY is not given, the correlations of the columns of datX are computed before testing.
#'In this case:
#'
#'- The differential correlation network can be expand from datX entites to a specific nodetype, by providing a value to \code{xTo}.
#'
#'If \code{xTo} or \code{yTo} or both is given, the columns of both datX and datY are required to use grinn ids for extended queries on the Grinn internal database, see \code{\link{convertToGrinnID}} for id conversion.
#'
#'If \code{xTo} = NULL and \code{yTo} = NULL, only the correlation network will be returned.
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Output includes correlation coefficients, pvalues and relation directions of each conditions, 
#'and the pvalues (pvalDiff) after testing. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references Langfelder P. and Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9:559 
#'@references Dudoit S., Yang YH., Callow MJ. and Speed TP. (2002) Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA, 12:111
#'@references Langfelder P. and Horvath S. Tutorials for the WGCNA package \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html}
#'@references Fukushima A. (2013) DiffCorr: an R package to analyze and visualize differential correlations in biological networks. Gene, 10;518(1):209-14.
#'@export
#'@seealso \code{\link{cor}}, \code{\link{corPvalueStudent}}, \code{\link{fetchGrinnNetwork}}, \pkg{\link{DiffCorr}}, \url{http://js.cytoscape.org/}
#'@examples
#'# Compute a differential correlation network of metabolites and expand to a grinn network of metabolite-protein
#'dummyX1 <- rbind(nodetype=rep("metabolite"),mtcars[1:16,])
#'colnames(dummyX1) <- c('G1.1','G27967','G371','G4.1',letters[1:7])
#'rownames(dummyX1)[-1] <- paste0(rep("normal_"),1:16)
#'dummyX2 <- rbind(nodetype=rep("metabolite"),mtcars[17:32,])
#'colnames(dummyX2) <- c('G1.1','G27967','G371','G4.1',letters[1:7])
#'rownames(dummyX2)[-1] <- paste0(rep("cancer_"),1:16)
#'result <- fetchDiffCorrGrinnNetwork(datX1=dummyX1, datX2=dummyX2, pDiff=0.05, method="spearman", returnAs="tab", yTo="protein")
#'library(igraph)
#'plot(graph.data.frame(result$edges[,1:2], directed=FALSE))
#'# Compute a differential correlation network of metabolites and proteins and expand to the grinn network of metabolite-pathway and protein-gene
#'dummyX1 <- rbind(nodetype=rep("metabolite"),mtcars[1:16,1:5])
#'colnames(dummyX1) <- c('G1.1','G27967','G371','G4.1','G16962')
#'rownames(dummyX1)[-1] <- paste0(rep("normal_"),1:16)
#'dummyX2 <- rbind(nodetype=rep("metabolite"),mtcars[17:32,1:5])
#'colnames(dummyX2) <- c('G1.1','G27967','G371','G4.1','G16962')
#'rownames(dummyX2)[-1] <- paste0(rep("cancer_"),1:16)
#'dummyY1 <- rbind(nodetype=rep("protein"),mtcars[1:16,6:10])
#'colnames(dummyY1) <- c('P28845','P08235','Q08AG9','P80365','P15538')
#'rownames(dummyY1)[-1] <- paste0(rep("normal_"),1:16)
#'dummyY2 <- rbind(nodetype=rep("protein"),mtcars[17:32,6:10])
#'colnames(dummyY2) <- c('P28845','P08235','Q08AG9','P80365','P15538')
#'rownames(dummyY2)[-1] <- paste0(rep("cancer_"),1:16)
#'result <- fetchDiffCorrGrinnNetwork(datX1=dummyX1, datX2=dummyX2, datY1=dummyY1, datY2=dummyY2, pDiff=0.05, method="spearman", returnAs="tab", xTo="pathway", yTo="gene")

fetchDiffCorrGrinnNetwork <- function(datX1=datX1,  datX2=datX2, datY1=NULL, datY2=NULL, pDiff=1e-4, method="spearman", returnAs="tab",
                                      xTo=NULL, yTo=NULL, filterSource=list()){
  corrnw = fetchDiffCorrNetwork(datX1=datX1,datX2=datX2,datY1=datY1,datY2=datY2,pDiff=pDiff,method=method,returnAs="tab")
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
        basicnw = fetchGrinnNetwork(txtInput=corrnw$nodes$id,from=nodetypes,
                                    to=xTo,filterSource=filterSource,dbXref="grinn")
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
      basicnw$edges$pvalDiff = 0
      basicnw$edges$condition = ""
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