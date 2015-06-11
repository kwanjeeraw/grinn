#'Compute a weighted correlation network and expand the network with information from grinn internal database
#'@description from input omics data e.g. normalized expression data or metabolomics data, it is a one step function to:
#'
#'1. Compute a weighted correlation network of input omics data using WGCNA functions \code{cor} and \code{corPvalueStudent}. 
#'The correlation coefficients are continuous values between -1 (negative correlation) and 1 (positive correlation), with numbers close to 1 or -1, meaning very closely correlated.
#'Then the correlation network is built by function \code{fetchCorrNetwork}.
#'
#'2. Expand the correlation network using information from grinn internal database.
#'The nodes of the correlation network are the keywords input to query the Grinn database.
#'Grinn internal database contains the networks of the following types that can get expanded to: 
#'metabolite-protein, metabolite-protein-gene, metabolite-pathway, protein-gene, protein-pathway and gene-pathway, see also \code{\link{fetchGrinnNetwork}}.
#'@usage fetchCorrGrinnNetwork(datNormX, datNormY, corrCoef, pval, method, returnAs, sourceTo, targetTo, filterSource)
#'@param datNormX data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities. 
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples e.g. normals, tumors. 
#'Require 'nodetype' at the first row to indicate the type of entities in each column. See below for details.
#'@param datNormY data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities.
#'Use the same format as \code{datNormX} or it can be NULL. See below for details.
#'@param corrCoef numerical value to define the minimum value of absolute correlation, from 0 to 1, to include edges in the output.
#'@param pval numerical value to define the maximum value of pvalues, to include edges in the output.
#'@param method string to define which correlation is to be used. It can be one of "pearson","kendall","spearman" (default), see \code{\link{cor}}.  
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape", default is "tab". "cytoscape" is the format used in Cytoscape.js
#'@param sourceTo string of node type. It can be one of "metabolite","protein","gene","pathway". See below for details.
#'@param targetTo string of node type. It can be one of "metabolite","protein","gene","pathway". By default, it will expand to pathways, see below for details.
#'@param filterSource string or list of pathway databases. The argument is required, if \code{sourceTo} or \code{targetTo = "pathway"}.
#'The argument value can be any of "SMPDB","KEGG","REACTOME" or combination of them e.g. list("KEGG","REACTOME").
#'@details datNormX and datNormY are matrices in which rows are samples and columns are entities.
#'If datNormY is given, then the correlations between the columns of datNormX and the columns of datNormY are computed.
#'In this case: 
#'
#'- The correlation network can be expand from datNormX entites to a specific node type, by providing a value to \code{sourceTo}
#'
#'- The correlation network can be expand from datNormY entites to a specific node type, by providing a value to \code{targetTo}
#'
#'Otherwise if datNormY is not given, the correlations of the columns of datNormX are computed. 
#'In this case:
#'
#'- The correlation network can be expand from datNormX entites to a specific node type, by providing a value to \code{targetTo} and leave \code{sourceTo = NULL}.
#'The column names of both datNormX and datNormY are required to use grinn ids. \code{convertToGrinnID} is provided for id conversion, see \code{\link{convertToGrinnID}}.
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references Langfelder P. and Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9:559 
#'@references Dudoit S., Yang YH., Callow MJ. and Speed TP. (2002) Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA, 12:111
#'@references Langfelder P. and Horvath S. Tutorials for the WGCNA package \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html}
#'@export
#'@seealso \code{\link{cor}}, \code{\link{corPvalueStudent}}, \code{\link{fetchCorrNetwork}}, \code{\link{fetchGrinnNetwork}}, \url{http://js.cytoscape.org/}
#'@examples
#'# Compute a correlation network of metabolites and expand to a grinn network of metabolite-protein
#'dummy <- rbind(nodetype=rep("metabolite"),t(mtcars))
#'colnames(dummy) <- c('X1.1','X27967','X371','X4.1',paste0('X',sample(400:22000, 28)))
#'result <- fetchCorrGrinnNetwork(datNormX=dummy, corrCoef=0.7, pval=1e-12, method="spearman", returnAs="tab", targetTo="protein")
#'library(igraph)
#'plot(graph.data.frame(result$edges[,1:2], directed=FALSE))
#'# Compute a correlation network of metabolites and proteins and expand to the grinn network of metabolite-pathway and protein-gene
#'dummyX <- rbind(nodetype=rep("metabolite"),t(mtcars)[,1:16])
#'colnames(dummyX) <- c('X1.1','X27967','X371','X4.1',paste0('X',sample(400:22000, 12)))
#'dummyY <- rbind(nodetype=rep("protein"),t(mtcars)[,17:32])
#'colnames(dummyY) <- c('P28845','P08235','Q08AG9','P80365',paste0('P',sample(10000:80000, 12)))
#'result <- fetchCorrGrinnNetwork(datNormX=dummyX, datNormY=dummyY, corrCoef=0.7, pval=1e-4, method="spearman", returnAs="json", sourceTo="pathway", targetTo="gene")

fetchCorrGrinnNetwork <- function(datNormX, datNormY=NULL, corrCoef=0.5, pval=1e-9, method="spearman", returnAs="tab", 
                                  sourceTo=NULL, targetTo=NULL, filterSource=list()){
  corrnw = fetchCorrNetwork(datNormX=datNormX,datNormY=datNormY,corrCoef=corrCoef,pval=pval,method=method,returnAs="tab")
  if(nrow(corrnw$nodes)>0){
    nodetypes = tolower(unique(corrnw$nodes$nodetype))
    if(length(nodetypes)>1){#if there are two data types
      #basicnw1 = fetchGrinnNetwork(txtInput=corrnw$nodes[which(tolower(corrnw$nodes$nodetype)==nodetypes[1]), 1],from=nodetypes[1],
      #                             to=nodetypes[2],filterSource=list(),dbXref="grinn") #relations between datasets
      #basicnw2 = fetchGrinnNetwork(txtInput=corrnw$nodes[which(tolower(corrnw$nodes$nodetype)==nodetypes[2]), 1],from=nodetypes[2],
      #                             to=nodetypes[1],filterSource=list(),dbXref="grinn") #relations between datasets
      #basicnw3 = rbind(basicnw1$edges,basicnw2$edges) #collect all edges
      #basicnw3 = basicnw3[duplicated(basicnw3[,1:2]),] #choose overlapping edges
      #basicnwNodes = rbind(basicnw1$nodes,basicnw2$nodes) #collect all nodes
      #basicnwNodes = basicnwNodes[duplicated(basicnwNodes[,1]),] #choose overlapping nodes
      if(!is.null(sourceTo)){
        basicnw4 = fetchGrinnNetwork(txtInput=corrnw$nodes[which(tolower(corrnw$nodes$nodetype)==nodetypes[1]), 1],from=nodetypes[1],
                                   to=sourceTo,filterSource=filterSource,dbXref="grinn") #relations sources to others
        basicnwNodes = basicnw4$nodes #collect all nodes
        #basicnwNodes = rbind(basicnwNodes,basicnw4$nodes) #collect all nodes
      }else{
        basicnw4 = data.frame()
        basicnwNodes = data.frame()
      }
      basicnw5 = fetchGrinnNetwork(txtInput=corrnw$nodes[which(tolower(corrnw$nodes$nodetype)==nodetypes[2]), 1],from=nodetypes[2],
                                   to=targetTo,filterSource=filterSource,dbXref="grinn") #relations from targets to others
      basicnwNodes = rbind(basicnwNodes,basicnw5$nodes) #collect all nodes
      basicnwEdges = rbind(basicnw4$edges,basicnw5$edges) #collect all edges
      #basicnwEdges = rbind(basicnw3,basicnw4$edges,basicnw5$edges) #collect all edges
      if(!is.null(basicnwEdges)){
        basicnwEdges = basicnwEdges[!duplicated(basicnwEdges[,1:ncol(basicnwEdges)]),] #remove duplicated edges
        basicnwNodes = basicnwNodes[!duplicated(basicnwNodes[,1]),] #remove duplicated nodes
        basicnw = list(nodes=basicnwNodes, edges=basicnwEdges)
      }else{
        basicnw = list(nodes=data.frame(),edges=data.frame())
      }
    }else{#if there is only one data type
      basicnw = fetchGrinnNetwork(txtInput=corrnw$nodes$id,from=nodetypes,to=targetTo,filterSource=filterSource,dbXref="grinn")
    }
    #collect node info
    corrattb = data.frame()
    for(i in 1:nrow(corrnw$nodes)){
      querystring <- paste0("MATCH (node:",corrnw$nodes$nodetype[i],") WHERE lower(node.GID) = lower('",corrnw$nodes$id[i],"') RETURN DISTINCT node")
      result <- curlRequestCypher(querystring)
      if(length(result)>0){
        corrattb = rbind(corrattb,data.frame(id=result[[1]]$data$GID,nodename=result[[1]]$data$name,nodetype=corrnw$nodes$nodetype[i],xref=paste0(unlist(result[[1]]$data$xref),collapse = "||")))
      }
    }
    if(nrow(basicnw$nodes)>0){
      cat("Formating and returning combined network ...\n")
      basicnw$edges$corr_coef = 1
      basicnw$edges$pval = 0
      basicnw$edges$direction = 0
      corrnw$edges$relsource = ""
      corrnw$nodes$xref = ""
      pair = rbind(basicnw$edges,corrnw$edges)
      if(nrow(corrattb)>0){attb = rbind(basicnw$nodes,corrattb,corrnw$nodes)}else{attb = rbind(basicnw$nodes,corrnw$nodes)}
      attb = attb[!duplicated(attb[,1]),]
      colnames(attb) = c("id","nodename","xref","nodetype")
      cat("Found ",nrow(pair)," relationships...\n")
    }else{#if only correlation network found
      cat("Formating and returning combined network ...\n")
      pair = corrnw$edges
      if(nrow(corrattb)>0){attb = rbind(corrattb,corrnw$nodes)}else{attb = corrnw$nodes}
      attb = attb[!duplicated(attb[,1]),]
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