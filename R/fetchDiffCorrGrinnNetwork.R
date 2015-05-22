#'Compute a differential correlation network and expand the network with information from grinn internal database
#'@description from input omics data e.g. normalized expression data or metabolomics data, it is a one step function to:
#'
#'1. Compute a differential correlation network of input omics data from two conditions, see \code{datNormX1}, \code{datNormX2}, \code{datNormY1}, \code{datNormY2}.
#'Correlation coefficients, pvalues and relation directions among entities in each condition are calculated using WGCNA functions \code{cor} and \code{corPvalueStudent}.
#'The correlation coefficients are continuous values between -1 (negative correlation) and 1 (positive correlation), with numbers close to 1 or -1, meaning very closely correlated.
#'Then correlation coefficients are test for differential correlations using Fisher's z-test based on \pkg{DiffCorr}.
#'The differential correlation network is created by function \code{fetchDiffCorrNetwork}.
#'
#'2. Expand the differential correlation network using information from grinn internal database.
#'The nodes of the differential correlation network are the keywords input to query the grinn database.
#'Grinn internal database contains the networks of the following types that can get expanded to: 
#'metabolite-protein, metabolite-protein-gene, metabolite-pathway, protein-gene, protein-pathway and gene-pathway, see also \code{\link{fetchGrinnNetwork}}.
#'@usage fetchDiffCorrGrinnNetwork(datNormX1,  datNormX2, datNormY1, datNormY2, pDiff, method, returnAs, sourceTo, targetTo, filterSource, organism)
#'@param datNormX1 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of one condition. 
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples. 
#'Require 'nodetype' at the first row to indicate the type of entities in each column. See below for details.
#'@param datNormX2 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of another condition. 
#'Use the same format as \code{datNormX1}.
#'@param datNormY1 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of one condition.
#'Use the same format as \code{datNormX1}. If there is only one type of dataset, it can be NULL. See below for details.
#'@param datNormY2 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of another condition.
#'Use the same format as \code{datNormX1}. If there is only one type of dataset, it can be NULL. See below for details.
#'@param pDiff numerical value to define the maximum value of pvalues (pvalDiff), to include edges in the output.
#'@param method string to define which correlation is to be used. It can be one of "pearson","kendall","spearman", see \code{\link{cor}}.  
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape", default is "tab". "cytoscape" is the format used in Cytoscape.js
#'@param sourceTo string of node type. It can be one of "metabolite","protein","gene","pathway". See below for details.
#'@param targetTo string of node type. It can be one of "metabolite","protein","gene","pathway". By default, it will expand to pathways, see below for details.
#'@param filterSource string or list of pathway databases. The argument is required, if \code{sourceTo} or \code{targetTo = "pathway"}.
#'The argument value can be any of "SMPDB","KEGG","REACTOME" or combination of them e.g. list("KEGG","REACTOME").
#'@param organism string of species in the following format: organism = "'species'". Default is "'Homo sapiens'".
#'@details To calculate the differential correlation network, require the input data from two conditions; 1 and 2. 
#'The input data are matrices in which rows are samples and columns are entities.
#'For each condition, if datNormY is given, then the correlations between the columns of datNormX and the columns of datNormY are computed before testing.
#'In this case: 
#'
#'- The differential correlation network can be expand from datNormX entites to a specific node type, by providing a value to \code{sourceTo}
#'
#'- The differential correlation network can be expand from datNormY entites to a specific node type, by providing a value to \code{targetTo}
#'
#'Otherwise if datNormY is not given, the correlations of the columns of datNormX are computed before testing.
#'In this case:
#'
#'- The correlation network can be expand from datNormX entites to a specific node type, by providing a value to \code{targetTo} and leave \code{sourceTo = NULL}.
#'The column names of the input data are required to use grinn ids. \code{convertToGrinnID} is provided for id conversion, see \code{\link{convertToGrinnID}}.
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Output includes correlation coefficients, pvalues and relation directions of each conditions, 
#'and the pvalues (pvalDiff) after testing. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references Langfelder P. and Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9:559 
#'@references Dudoit S., Yang YH., Callow MJ. and Speed TP. (2002) Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA, 12:111
#'@references Langfelder P. and Horvath S. Tutorials for the WGCNA package \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html}
#'@references Fukushima A. (2013) DiffCorr: an R package to analyze and visualize differential correlations in biological networks. Gene, 10;518(1):209-14.
#'@export
#'@seealso \code{\link{cor}}, \code{\link{corPvalueStudent}}, \code{\link{fetchDiffCorrNetwork}}, \code{\link{fetchGrinnNetwork}}, \pkg{\link{DiffCorr}}, \url{http://js.cytoscape.org/}
#'@examples
#'# Compute a differential correlation network of metabolites and expand to a grinn network of metabolite-protein
#'dummyX1 <- rbind(nodetype=rep("metabolite"),mtcars[1:16,])
#'colnames(dummyX1) <- c('X1.1','X27967','X371','X4.1',letters[1:7])
#'rownames(dummyX1)[-1] <- paste0(rep("normal_"),1:16)
#'dummyX2 <- rbind(nodetype=rep("metabolite"),mtcars[17:32,])
#'colnames(dummyX2) <- c('X1.1','X27967','X371','X4.1',letters[1:7])
#'rownames(dummyX2)[-1] <- paste0(rep("cancer_"),1:16)
#'result <- fetchDiffCorrGrinnNetwork(datNormX1=dummyX1, datNormX2=dummyX2, pDiff=0.05, method="spearman", returnAs="tab", targetTo="protein")
#'library(igraph)
#'plot(graph.data.frame(result$edges[,1:2], directed=FALSE))
#'# Compute a differential correlation network of metabolites and proteins and expand to the grinn network of metabolite-pathway and protein-gene
#'dummyX1 <- rbind(nodetype=rep("metabolite"),mtcars[1:16,1:5])
#'colnames(dummyX1) <- c('X1.1','X27967','X371','X4.1','X16962')
#'rownames(dummyX1)[-1] <- paste0(rep("normal_"),1:16)
#'dummyX2 <- rbind(nodetype=rep("metabolite"),mtcars[17:32,1:5])
#'colnames(dummyX2) <- c('X1.1','X27967','X371','X4.1','X16962')
#'rownames(dummyX2)[-1] <- paste0(rep("cancer_"),1:16)
#'dummyY1 <- rbind(nodetype=rep("protein"),mtcars[1:16,6:10])
#'colnames(dummyY1) <- c('P28845','P08235','Q08AG9','P80365','P15538')
#'rownames(dummyY1)[-1] <- paste0(rep("normal_"),1:16)
#'dummyY2 <- rbind(nodetype=rep("protein"),mtcars[17:32,6:10])
#'colnames(dummyY2) <- c('P28845','P08235','Q08AG9','P80365','P15538')
#'rownames(dummyY2)[-1] <- paste0(rep("cancer_"),1:16)
#'result <- fetchDiffCorrGrinnNetwork(datNormX1=dummyX1, datNormX2=dummyX2, datNormY1=dummyY1, datNormY2=dummyY2, pDiff=0.05, method="spearman", returnAs="tab", sourceTo="pathway", targetTo="gene")

fetchDiffCorrGrinnNetwork <- function(datNormX1=datNormX1,  datNormX2=datNormX2, datNormY1=NULL, datNormY2=NULL, pDiff=1e-4, method="spearman", returnAs="tab",
                                      sourceTo=NULL, targetTo=NULL, filterSource=list(), organism="'Homo sapiens'"){
  corrnw = fetchDiffCorrNetwork(datNormX1=datNormX1,datNormX2=datNormX2,datNormY1=datNormY1,datNormY2=datNormY2,pDiff=pDiff,method=method,returnAs="tab")
  if(nrow(corrnw$nodes)>0){
    nodetypes = tolower(unique(corrnw$nodes$nodetype))
    if(length(nodetypes)>1){#if there are two data types
      basicnw1 = fetchGrinnNetwork(txtInput=corrnw$nodes[which(tolower(corrnw$nodes$nodetype)==nodetypes[1]), 1],from=nodetypes[1],
                                   to=nodetypes[2],filterSource=list(),dbXref="grinn") #relations between datasets
      basicnw2 = fetchGrinnNetwork(txtInput=corrnw$nodes[which(tolower(corrnw$nodes$nodetype)==nodetypes[2]), 1],from=nodetypes[2],
                                   to=nodetypes[1],filterSource=list(),dbXref="grinn") #relations between datasets
      basicnw3 = rbind(basicnw1$edges,basicnw2$edges) #collect all edges
      basicnw3 = basicnw3[duplicated(basicnw3[,1:2]),] #choose overlapping edges
      basicnwNodes = rbind(basicnw1$nodes,basicnw2$nodes) #collect all nodes
      basicnwNodes = basicnwNodes[duplicated(basicnwNodes[,1]),] #choose overlapping nodes
      if(!is.null(sourceTo)){
        basicnw4 = fetchGrinnNetwork(txtInput=corrnw$nodes[which(tolower(corrnw$nodes$nodetype)==nodetypes[1]), 1],from=nodetypes[1],
                                    to=sourceTo,filterSource=filterSource,dbXref="grinn") #relations sources to others
        basicnwNodes = rbind(basicnwNodes,basicnw4$nodes) #collect all nodes
      }else{
        basicnw4 = data.frame()
      }
      basicnw5 = fetchGrinnNetwork(txtInput=corrnw$nodes[which(tolower(corrnw$nodes$nodetype)==nodetypes[2]), 1],from=nodetypes[2],
                                   to=targetTo,filterSource=filterSource,dbXref="grinn") #relations from targets to others
      basicnwNodes = rbind(basicnwNodes,basicnw5$nodes) #collect all nodes
      
      basicnwEdges = rbind(basicnw3,basicnw4$edges,basicnw5$edges) #collect all edges
      if(!is.null(basicnwEdges)){
        basicnwEdges = basicnwEdges[!duplicated(basicnwEdges[,1:ncol(basicnwEdges)]),] #remove duplicated edges
        basicnwNodes = basicnwNodes[!duplicated(basicnwNodes[,1]),] #remove duplicated nodes
        basicnw = list(nodes=basicnwNodes, edges=basicnwEdges)
      }else{
        basicnw = list(nodes=data.frame(),edges=data.frame())
      }
    }else{#if there is only one data type
      basicnw = fetchGrinnNetwork(txtInput=corrnw$nodes$id,from=nodetypes,
                                  to=targetTo,filterSource=filterSource,dbXref="grinn")
    }
    if(nrow(basicnw$nodes)>0){
      cat("Formating and returning combined network ...\n")
      basicnw$edges$corr_coef = 1
      basicnw$edges$pval = 0
      basicnw$edges$direction = 0
      basicnw$edges$pvalDiff = 0
      basicnw$edges$condition = ""
      corrnw$edges$relsource = ""
      corrnw$nodes$xref = ""
      pair = rbind(basicnw$edges,corrnw$edges)
      attb = rbind(basicnw$nodes,corrnw$nodes)
      attb = attb[!duplicated(attb[,1]),]
      colnames(attb) = c("id","nodename","xref","nodetype")
      cat("Found ",nrow(pair)," relationships...\n")
    }else{#if only correlation network found
      cat("Formating and returning combined network ...\n")
      pair = corrnw$edges
      attb = corrnw$nodes
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