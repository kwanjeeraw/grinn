#'Compute a weighted correlation network
#'@description from the input data e.g. normalized expression data, build the correlation network based on \pkg{WGCNA}. 
#'Correlation coefficients, pvalues and relation directions are calculated using WGCNA functions \code{cor} and \code{corPvalueStudent}.
#'The correlation coefficients are continuous values between -1 (negative correlation) and 1 (positive correlation), with numbers close to 1 or -1, meaning very closely correlated.
#'datNormX and datNormY are matrices in which rows are samples and columns are entities.
#'@usage fetchCorrNetwork(datNormX, datNormY, corrCoef, pval, method, returnAs)
#'@param datNormX data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities. 
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples e.g. normals, tumors. 
#'Require 'nodetype' at the first row to indicate the type of entities in each column. See below for details.
#'@param datNormY data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities.
#'Use the same format as \code{datNormX}. If there is only one type of dataset, \code{datNormY} must be \code{datNormY = NULL}. See below for details.
#'@param corrCoef numerical value to define the minimum value of absolute correlation, from 0 to 1, to include edges in the output.
#'@param pval numerical value to define the maximum value of pvalues, to include edges in the output.
#'@param method string to define which correlation is to be used. It can be one of "pearson","kendall","spearman", see \code{\link{cor}}.  
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape". "cytoscape" is the format used in Cytoscape.js
#'@details datNormX and datNormY are matrices in which rows are samples and columns are entities.
#'If datNormY is given, then the correlations between the columns of datNormX and the columns of datNormY are computed.
#'Otherwise if datNormY is not given, the correlations of the columns of datNormX are computed. 
#'If grinn functions will be used in further analyses, the column names of both datNormX and datNormY are suggested to use grinn ids. 
#'\code{convertToGrinnID} is provided for id conversion, see \code{\link{convertToGrinnID}}.
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references Langfelder P. and Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9:559 
#'@references Dudoit S., Yang YH., Callow MJ. and Speed TP. (2002) Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA, 12:111
#'@references Langfelder P. and Horvath S. Tutorials for the WGCNA package \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html}
#'@export
#'@seealso \code{\link{cor}}, \code{\link{corPvalueStudent}}, \url{http://js.cytoscape.org/}
#'@examples
#'# Compute a correlation network of metabolites
#'dummy <- rbind(nodetype=rep("metabolite"),t(mtcars))
#'colnames(dummy) <- c('X1.1','X27967','X371','X4.1',paste0('X',sample(400:22000, 28)))
#'result <- fetchCorrNetwork(datNormX=dummy, datNormY=NULL, corrCoef=0.7, pval=1e-12, method="spearman", returnAs="tab")
#'library(igraph)
#'plot(graph.data.frame(result$edges[,1:2], directed=FALSE))
#'# Compute a correlation network of metabolites and proteins
#'dummyX <- rbind(nodetype=rep("metabolite"),t(mtcars)[,1:16])
#'colnames(dummyX) <- c('X1.1','X27967','X371','X4.1',paste0('X',sample(400:22000, 12)))
#'dummyY <- rbind(nodetype=rep("protein"),t(mtcars)[,17:32])
#'colnames(dummyY) <- c('P28845','P08235','Q08AG9','P80365',paste0('P',sample(10000:80000, 12)))
#'result <- fetchCorrNetwork(datNormX=dummyX, datNormY=dummyY, corrCoef=0.7, pval=1e-4, method="spearman", returnAs="tab")

fetchCorrNetwork <- function(datNormX, datNormY, corrCoef, pval, method, returnAs){
  #calculate correlation network adjacency from normalized data
  if(!("nodetype" %in% rownames(datNormX))){
    stop("can't define type of entities, missing nodetype")
  }

  #create data frame from matrix
  mtTodf <- function(mt) {
    upmt <- upper.tri(mt)
    data.frame(from = rownames(mt)[row(mt)[upmt]], to = rownames(mt)[col(mt)[upmt]], val = t(mt)[upmt], stringsAsFactors = F)
  }
  #format reltype
  maptype = function(x){
    rtype = paste0(Hmisc::capitalize(gsub(".+TYPE_","",x[1])),"_",Hmisc::capitalize(gsub(".+TYPE_","",x[2])))
  }
  corAdj = getCorrAdjacency(datNormX,datNormY,method)
  #filter by corr_coef and pval
  ix = which(abs(corAdj$corr_coef) > corrCoef)
  iy = which(corAdj$pval < pval)
  it = intersect(ix, iy)
  if(length(it)>0){
    cat("Formating and returning correlation network ...\n")
    corDF = corAdj[it,]
    corDF$reltype = apply(corDF,1,maptype)    
    lnodes = unique(c(as.character(corDF$source),as.character(corDF$target)))
    attb = data.frame(id=gsub("_TYPE.+","",lnodes),name=gsub("_TYPE.+","",lnodes),nodetype=Hmisc::capitalize(gsub(".+TYPE_","",lnodes)), stringsAsFactors = F)
    if(length(it) != 1){ 
      st = apply(corDF[,1:2],2,function(x) gsub("_TYPE.+","",x)) 
    }else{
      st = t(as.matrix(apply(corDF[,1:2],2,function(x) gsub("_TYPE.+","",x))))
    }
    pair = data.frame(st,corDF[,3:ncol(corDF)], stringsAsFactors = F)
    pair$relname = "CORRELATION"
    colnames(attb) = c("id","nodename","nodetype")
    cat("Found ",nrow(pair)," correlations...\n")
  }else{# if no network found
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