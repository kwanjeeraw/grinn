#'Compute a differential correlation network
#'@description take the input data from two conditions e.g. normalized gene expression data from normal and cancer cells, 
#'and then compute the differential correlation network based on \pkg{DiffCorr}. 
#'Correlation coefficients, pvalues and relation directions among entities in each condition are calculated using WGCNA functions \code{cor} and \code{corPvalueStudent}.
#'The correlation coefficients are continuous values between -1 (negative correlation) and 1 (positive correlation), with numbers close to 1 or -1, meaning very closely correlated.
#'Then correlation coefficients are test for differential correlations using Fisher's z-test based on \pkg{DiffCorr}.
#'@usage fetchDiffCorrNetwork(datNormX1, datNormX2, datNormY1, datNormY2, pDiff, method, returnAs)
#'@param datNormX1 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of one condition. 
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples. 
#'Require 'nodetype' at the first row to indicate the type of entities in each column. See below for details.
#'@param datNormX2 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of another condition. 
#'Use the same format as \code{datNormX1}.
#'@param datNormY1 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of one condition.
#'Use the same format as \code{datNormX1}. If there is only one type of dataset, \code{datNormY1} must be \code{datNormY1 = NULL}. See below for details.
#'@param datNormY2 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of another condition.
#'Use the same format as \code{datNormX1}. If there is only one type of dataset, \code{datNormY2} must be \code{datNormY2 = NULL}. See below for details.
#'@param pDiff numerical value to define the maximum value of pvalues (pvalDiff), to include edges in the output.
#'@param method string to define which correlation is to be used. It can be one of "pearson","kendall","spearman", see \code{\link{cor}}.  
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape". "cytoscape" is the format used in Cytoscape.js
#'@details To calculate the differential correlation network, require the input data from two conditions; 1 and 2. The input data are matrices in which rows are samples and columns are entities.
#'For each condition, if datNormY is given, then the correlations between the columns of datNormX and the columns of datNormY are computed.
#'Otherwise if datNormY is not given, the correlations of the columns of datNormX are computed. Then correlation coefficients are test for significant correlation pairs.
#'If grinn functions will be used in further analyses, the column names of both datNormX and datNormY are suggested to use grinn ids.
#'\code{convertToGrinnID} is provided for id conversion, see \code{link{convertToGrinnID}}.
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Output includes correlation coefficients, pvalues and relation directions of each conditions, 
#'and the pvalues (pvalDiff) after testing. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references Langfelder P. and Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9:559 
#'@references Dudoit S., Yang YH., Callow MJ. and Speed TP. (2002) Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA, 12:111
#'@references Langfelder P. and Horvath S. Tutorials for the WGCNA package \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html}
#'@references Fukushima A. (2013) DiffCorr: an R package to analyze and visualize differential correlations in biological networks. Gene, 10;518(1):209-14.
#'@export
#'@seealso \code{\link{cor}}, \code{\link{corPvalueStudent}}, \pkg{\link{DiffCorr}}, \url{http://js.cytoscape.org/}
#'@examples
#'# Compute a differential correlation network of metabolites
#'dummyX1 <- rbind(nodetype=rep("metabolite"),mtcars[1:16,])
#'colnames(dummyX1) <- c('X1.1','X27967','X371','X4.1',letters[1:7])
#'rownames(dummyX1)[-1] <- paste0(rep("normal_"),1:16)
#'dummyX2 <- rbind(nodetype=rep("metabolite"),mtcars[17:32,])
#'colnames(dummyX2) <- c('X1.1','X27967','X371','X4.1',letters[1:7])
#'rownames(dummyX2)[-1] <- paste0(rep("cancer_"),1:16)
#'result <- fetchDiffCorrNetwork(datNormX1=dummyX1, datNormX2=dummyX2, datNormY1=NULL, datNormY2=NULL, pDiff=0.05, method="spearman", returnAs="tab")
#'library(igraph)
#'plot(graph.data.frame(result$edges[,1:2], directed=FALSE))
#'# Compute a differential correlation network of metabolites and proteins
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
#'result <- fetchDiffCorrNetwork(datNormX1=dummyX1, datNormX2=dummyX2, datNormY1=dummyY1, datNormY2=dummyY2, pDiff=0.05, method="spearman", returnAs="tab")

fetchDiffCorrNetwork <- function(datNormX1, datNormX2, datNormY1, datNormY2, pDiff, method, returnAs){
  #format reltype
  maptype = function(x){
    rtype = paste0(Hmisc::capitalize(gsub(".+TYPE_","",x[1])),"_",Hmisc::capitalize(gsub(".+TYPE_","",x[2])))
  }
  n1 = nrow(datNormX1) - 1 #number of samples in condition 1, remove first row = it is nodetype
  n2 = nrow(datNormX2) - 1 #number of samples in condition 2, remove first row = it is nodetype
  corrnw1 = getCorrAdjacency(datNormX=datNormX1,datNormY=datNormY1,method=method)
  corrnw2 = getCorrAdjacency(datNormX=datNormX2,datNormY=datNormY2,method=method)
  #transform corr_coef for testing, using equation from DiffCorr
  z1 = atanh(corrnw1$corr_coef)
  z2 = atanh(corrnw2$corr_coef)
  dz = (z1 - z2)/sqrt(1/(n1 - 3) + (1/(n2 - 3)))
  difp = 2 * (1 - pnorm(abs(dz)))
  corrnw1$pvalDiff = difp
  corrnw1$condition = "CONDITION1"
  corrnw2$pvalDiff = difp
  corrnw2$condition = "CONDITION2"
  ind = which(difp < pDiff)  
  if(length(ind)>0){
    cat("Formating and returning differential correlation network ...\n")
    corDF = corrnw1[ind,]
    corDF$reltype = apply(corDF,1,maptype)    
    lnodes = unique(c(as.character(corDF$source),as.character(corDF$target)))
    attb = data.frame(id=gsub("_TYPE.+","",lnodes),name=gsub("_TYPE.+","",lnodes),nodetype=Hmisc::capitalize(gsub(".+TYPE_","",lnodes)), stringsAsFactors = F)
    if(length(ind) != 1){ 
      st = apply(corDF[,1:2],2,function(x) gsub("_TYPE.+","",x)) 
    }else{
      st = t(as.matrix(apply(corDF[,1:2],2,function(x) gsub("_TYPE.+","",x))))
    }
    pair = data.frame(st,corDF[,3:ncol(corDF)], stringsAsFactors = F)
    tmpPair = data.frame(pair[,1:2],corrnw2[ind,3:ncol(corrnw2)],reltype=pair[,"reltype"])
    pair = rbind(pair,tmpPair)
    pair$relname = "DIFF_CORRELATION"
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