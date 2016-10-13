#'Compute partial correlation network
#'@description compute the partial correlation network of entities from raw or quantified data, see details.
#'@usage fetchPtCorrNetwork(datX,  corrCoef, pval, alpha, epsilon, matrix.completion, returnAs)
#'@param datX data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities. 
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples e.g. normals, tumors. 
#'Require 'nodetype' at the first row to indicate the type of entities in each column. See below for details.
#'@param corrCoef numerical value to define the minimum value of absolute correlation, from 0 to 1, to include edges in the output.
#'@param pval numerical value to define the maximum value of pvalues, to include edges in the output.
#'@param alpha a numeric value specifying significance level of each test used in \code{\link{qpAvgNrr}}.
#'@param epsilon a numeric value specifying the maximum cutoff value of the non-rejection rate met by the edges that are included in the qp-graph, see \code{\link{qpGraph}}.
#'@param matrix.completion a string specifying algorithm to employ in the matrix completion operations used in \code{\link{qpPAC}}
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape". "cytoscape" is the format used in Cytoscape.js
#'@details The function wraps around the functions of \pkg{\link{qpgraph}}.
#'Partial correlation coefficients, p-values and correlation directions are calculated.
#'The partial correlation coefficients are continuous values between -1 (negative correlation) and 1 (positive correlation), with numbers close to 1 or -1, meaning very closely correlated.
#'The datX is matrix in which rows are samples and columns are entities.
#'If grinn functions will be used in further analyses, the column names of datX are suggested to use grinn ids. 
#'\code{convertToGrinnID} is provided for id conversion, see \code{\link{convertToGrinnID}}.
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@seealso \pkg{\link{qpgraph}}, \code{\link{qpAvgNrr}}, \code{\link{qpGraph}}, \code{\link{qpPAC}}
#'@references Castelo R. and Roverato A. (2006) A robust procedure for Gaussian graphical model search from microarray data with p larger than n. J. Mach. Learn. Res., 7:2621-2650.
#'@references Castelo R. and Roverato A. (2009) Reverse engineering molecular regulatory networks from microarray data with qp-graphs. J Comput Biol, 16(2), pp. 213-27.
#'@examples
#'#dt = read.csv("~/Documents/grinn_sample/lung_miyamoto/metAdj.txt",sep="\t",header=TRUE,row.names=1)
#'#nw = fetchPtCorrNetwork(datX=dt, corrCoef=0.5, pval=0.05, returnAs="tab")
#'# Compute a par-correlation network of metabolites
#'#dummy <- rbind(nodetype=rep("metabolite"),t(mtcars))
#'#colnames(dummy) <- c('G1.1','G27967','G371','G4.1',paste0('G',sample(400:22000, 28)))
#'#result <- fetchPtCorrNetwork(datX=dummy, corrCoef=0.7, pval=0.05, returnAs="tab")
#'#library(igraph)
#'#plot(graph.data.frame(result$edges[,1:2], directed=FALSE))
fetchPtCorrNetwork <- function (datX,  corrCoef, pval, alpha, epsilon, matrix.completion, returnAs) 
 {
  tmparg <- try(matrix.completion <- match.arg(toupper(matrix.completion), c("IPF","HTF"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'matrix.completion' is not valid, choose one from the list: IPF,HTF")
  }
  #calculate correlation network adjacency from normalized data
  if(!("nodetype" %in% rownames(datX))){
    stop("can't define type of entities, missing nodetype")
  }
  cat("Formating row.names of input data frame ...\n")
  nodetype = datX[1,]
  colnames(datX) = paste0(colnames(datX),"_TYPE_",nodetype)
  datX = apply(datX[-1,],2,as.numeric)
  datX = t(datX)
  gc()
  cat("Computing partial correlation ...\n")
  nrr.estimates = qpgraph::qpAvgNrr(datX, alpha=alpha)
  g = qpgraph::qpGraph(nrr.estimates, epsilon=epsilon)
  pac.estimates = qpgraph::qpPAC(datX, g=g@g, matrix.completion=matrix.completion)
  #format output
  nRow = nrow(pac.estimates$R)
  nNames = dimnames(pac.estimates$R)[[1]]
  rowMat = matrix(c(1:nRow), nRow, nRow, byrow = TRUE)
  colMat = matrix(c(1:nRow), nRow, nRow)
  dstRows = as.dist(rowMat)
  dstCols = as.dist(colMat)
  network = data.frame(source = as.character(nNames[dstRows]), target = as.character(nNames[dstCols]), corr_coef = pac.estimates$R[lower.tri(pac.estimates$R)],
                       pval = pac.estimates$P[lower.tri(pac.estimates$P)], direction = sign(pac.estimates$R[lower.tri(pac.estimates$R)]), relname="PAR_CORRELATION", stringsAsFactors = FALSE)
  network = network[!is.na(network$pval),]
  network = network[abs(network$corr_coef) > corrCoef, ]
  network = network[network$pval < pval, ]
  if(nrow(network)>0){
    cat("Formating and returning network ...\n")
    lnodes = unique(c(as.character(network$source),as.character(network$target)))
    attb = data.frame(id=gsub("_TYPE.+","",lnodes),nodename=gsub("_TYPE.+","",lnodes),nodetype=Hmisc::capitalize(gsub(".+TYPE_","",lnodes)), stringsAsFactors = FALSE)
    st = apply(network[,1:2],2,function(x) gsub("_TYPE.+","",x)) 
    pair = data.frame(st,network[,3:ncol(network)], stringsAsFactors = FALSE)
    cat("Found ",nrow(attb)," nodes and ",nrow(pair)," edges...\n")
  }else{# if no network found
    print("Returning no data...")
    pair = data.frame()
    attb = data.frame()
  }
  out = switch(returnAs,
               tab = list(nodes=attb, edges=pair),
               json = list(nodes=jsonlite::toJSON(attb), edges=jsonlite::toJSON(pair)),
               stop("incorrect return type"))
}