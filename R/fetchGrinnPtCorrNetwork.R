#'Combine a grinn network queried from Grinn internal database to a partial correlation network
#'@description  from the list of keywords and input omics data e.g. normalized expression data or metabolomics data, it is a one step function to:
#'
#'1. Build an integrated network (grinn network) by connecting these keywords to a specified node type, see \code{\link{fetchGrinnNetwork}}.
#'The keywords can be any of these node types: metabolite, protein, gene and pathway.
#'The Grinn internal database contains the networks of the following types that can be quried: 
#'metabolite-protein, metabolite-protein-gene, metabolite-pathway, protein-gene, protein-pathway and gene-pathway. 
#'
#'2. Compute a partial correlation network of input omics data, see \code{datX}.
#'Correlation coefficients, pvalues and relation directions are calculated using qpgraph functions.
#'The correlation coefficients are continuous values between -1 (negative correlation) and 1 (positive correlation), with numbers close to 1 or -1, meaning very closely correlated.
#'
#'3. Combine the grinn network to the correlation network.
#'@usage fetchGrinnCorrNetwork(txtInput, from, to, filterSource, returnAs, dbXref, datX, corrCoef, pval, alpha, epsilon, matrix.completion)
#'@param txtInput list of keywords containing keyword ids e.g. txtInput = list('id1', 'id2'). 
#'The keyword ids are from the specified database, see \code{dbXref}. Default is grinn id e.g. G371.
#'@param from string of start node. It can be one of "metabolite","protein","gene","pathway".
#'@param to string of end node. It can be one of "metabolite","protein","gene","pathway".
#'@param filterSource string or list of pathway databases. The argument is required, if \code{from} or \code{to = "pathway"}, see \code{from} and \code{to}.
#'The argument value can be any of "SMPDB","KEGG","REACTOME" or combination of them e.g. list("KEGG","REACTOME").  
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape", default is "tab". "cytoscape" is the format used in Cytoscape.js
#'@param dbXref string of database name. Specify the database name used for the txtInput ids, see \code{txtInput}. 
#'It can be one of "grinn","chebi","kegg","pubchem","inchi","hmdb","smpdb","reactome","uniprot","ensembl","entrezgene". Default is "grinn".
#'If pubchem is used, it has to be pubchem SID (substance ID).
#'@param datX data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities. 
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples e.g. normals, tumors. 
#'Require 'nodetype' at the first row to indicate the type of entities in each column. See below for details.
#'@param corrCoef numerical value to define the minimum value of absolute correlation, from 0 to 1, to include edges in the output.
#'@param pval numerical value to define the maximum value of pvalues, to include edges in the output.
#'@param alpha a numeric value specifying significance level of each test used in \code{\link{qpAvgNrr}}.
#'@param epsilon a numeric value specifying the maximum cutoff value of the non-rejection rate met by the edges that are included in the qp-graph, see \code{\link{qpGraph}}.
#'@param matrix.completion a string specifying algorithm to employ in the matrix completion operations used in \code{\link{qpPAC}}
#'@details datX is matrix in which rows are samples and columns are entities. The correlations of the columns of datX are computed. 
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references Castelo R. and Roverato A. (2006) A robust procedure for Gaussian graphical model search from microarray data with p larger than n. J. Mach. Learn. Res., 7:2621-2650.
#'@references Castelo R. and Roverato A. (2009) Reverse engineering molecular regulatory networks from microarray data with qp-graphs. J Comput Biol, 16(2), pp. 213-27.
#'@export
#'@seealso \pkg{\link{qpgraph}}, \code{\link{qpAvgNrr}}, \code{\link{qpGraph}}, \code{\link{qpPAC}}, \code{\link{fetchGrinnNetwork}}, \url{http://js.cytoscape.org/}
#'@examples
#'# Create metabolite-gene network from the list of metabolites using grinn ids and combine the grinn network to a correlation network of metabolites
#'kw <- c('G160','G300','G371')
#'dummy <- rbind(nodetype=rep("metabolite"),t(mtcars))
#'colnames(dummy) <- c('G1.1','G27967','G371','G4.1',paste0('G',sample(400:22000, 28)))
#'result <- fetchGrinnPtCorrNetwork(txtInput=kw, from="metabolite", to="gene", datX=dummy, corrCoef=0.7, pval=0.05)
#'library(igraph)
#'plot(graph.data.frame(result$edges[,1:2], directed=FALSE))

fetchGrinnPtCorrNetwork <- function(txtInput, from, to, filterSource=list(), returnAs="tab", dbXref="grinn", datX, corrCoef=0.5, pval=0.05, alpha=0.05, epsilon=0.5, matrix.completion="IPF"){
  basicnw = fetchGrinnNetwork(txtInput=txtInput,from=from,to=to,filterSource=filterSource,dbXref=dbXref)
  corrnw = fetchPtCorrNetwork(datX=datX,corrCoef=corrCoef,pval=pval,alpha=alpha, epsilon=epsilon, matrix.completion=matrix.completion,returnAs="tab")
  if(nrow(corrnw$nodes)>0){
    #collect node info
    corrattb = data.frame()
    corrattb = plyr::ldply (apply(corrnw$nodes, MARGIN = 1, FUN=getNodeInfo, x = "id", y = "nodetype")) #format nodelist
    corrnw$edges$source = lapply(corrnw$edges$source, FUN=formatId, y = corrattb) #format edgelist
    corrnw$edges$target = lapply(corrnw$edges$target, FUN=formatId, y = corrattb) #format edgelist
  }
  if(nrow(basicnw$nodes)>0 && nrow(corrnw$nodes)>0){
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
  }else if(nrow(basicnw$nodes)>0 && nrow(corrnw$nodes)==0){
    cat("Formating and returning combined network ...\n")
    pair = basicnw$edges
    attb = basicnw$nodes
    cat("Found ",nrow(pair)," relationships...\n")
  }else if(nrow(basicnw$nodes)==0 && nrow(corrnw$nodes)>0){
    cat("Formating and returning combined network ...\n")
    pair = corrnw$edges
    corrnw$nodes$xref = ""
    corrnw$nodes$gid = corrnw$nodes$id #same ids
    if(nrow(corrattb)>0){attb = rbind(corrattb,corrnw$nodes)}else{attb = corrnw$nodes}
    attb = attb[!duplicated(attb[,2]),]
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