#' \code{fetchGrinnDiffCorrNetwork} combine a grinn network queried from Grinn internal database to a differential correlation network
#'@description  from the list of keywords and input omics data e.g. normalized expression data or metabolomics data, it is a one step function to:
#'1. Build an integrated network (grinn network) by connecting these keywords to a specified node type, see \code{\link{fetchGrinnNetwork}}.
#'The keywords can be any of these node types: metabolite, protein, gene and pathway.
#'Grinn internal database contains the networks of the following types that can be quried: 
#'metabolite-protein, metabolite-protein-gene, metabolite-pathway, protein-gene, protein-pathway and gene-pathway. 
#'2. Compute a differential correlation network of input omics data from two conditions, see \code{datNormX1}, \code{datNormX2}, \code{datNormY1}, \code{datNormY2}.
#'Correlation coefficients, pvalues and relation directions among entities in each condition are calculated using \code{WGCNA::cor} and \code{WGCNA::corPvalueStudent}.
#'The correlation coefficients are continuous values between -1 (negative correlation) and 1 (positive correlation), with numbers close to 1 or -1, meaning very closely correlated.
#'Then correlation coefficients are test for differential correlations using Fisher's z-test based on \pkg{DiffCorr}.
#'The differential correlation network is created by function \code{fetchDiffCorrNetwork}.
#'3. Combine the grinn network to the correlation network.
#'@usage fetchGrinnDiffCorrNetwork(txtInput, from, to, filterSource, returnAs, dbXref, organism, datNormX1, datNormX2, datNormY1, datNormY2, pDiff, method)
#'@param txtInput list of keywords containing keyword ids e.g. txtInput = list('id1', 'id2'). 
#'The keyword ids are from the specified database, see \code{dbXref}. Default is grinn id e.g. X371.
#'@param from string of start node. It can be one of "metabolite","protein","gene","pathway".
#'@param to string of end node. It can be one of "metabolite","protein","gene","pathway".
#'@param filterSource string or list of pathway databases. The argument is required, if \code{from} or \code{to} = "pathway", see \code{from} and \code{to}.
#'The argument value can be any of "SMPDB","KEGG","REACTOME" or combination of them e.g. list("KEGG","REACTOME").  
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape", default is "tab". "cytoscape" is the format used in Cytoscape.js
#'@param dbXref string of database name. Specify the database name used for the txtInput ids, see \code{txtInput}. 
#'It can be one of "grinn","chebi","kegg","pubchem","inchi","hmdb","smpdb","reactome","uniprot","ensembl","entrezgene". Default is "grinn".
#'If pubchem is used, it has to be pubchem SID (substance ID).
#'@param organism string of species in the following format: organism = "'species'". Default is "'Homo sapiens'".
#'@param datNormX1 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of one condition. 
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples e.g. normals, tumors. 
#'Require 'nodetype' at the first row to indicate the type of entities in each column. See below for details.
#'@param datNormX2 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of another condition. 
#'Use the same format as \code{datNormX1}.
#'@param datNormY1 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of one condition.
#'Use the same format as \code{datNormX1}. If there is only one type of dataset, it can be NULL. See below for details.
#'@param datNormY2 data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities of another condition.
#'Use the same format as \code{datNormX1}. If there is only one type of dataset, it can be NULL. See below for details.
#'@param pDiff numerical value to define the maximum value of pvalues (pvalDiff), to include edges in the output.
#'@param method string to define which correlation is to be used. It can be one of "pearson","kendall","spearman", see \code{WGCNA::cor}.  
#'@details To calculate the differential correlation network, require the input data from two conditions; 1 and 2. 
#'The input data are matrices in which rows are samples and columns are entities.
#'For each condition, if datNormY is given, then the correlations between the columns of datNormX and the columns of datNormY are computed.
#'Otherwise if datNormY is not given, the correlations of the columns of datNormX are computed. Then correlation coefficients are test for significant correlation pairs.
#'The column names of the input data are required to use grinn ids. \code{convertToGrinnID} is provided for id conversion, see \code{convertToGrinnID}.
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Output includes correlation coefficients, pvalues and relation directions of each conditions, 
#'and the pvalues (pvalDiff) after testing. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references 
#'Langfelder P. and Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9:559 
#'Dudoit S., Yang YH., Callow MJ. and Speed TP. (2002) Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA, 12:111
#'Langfelder P. and Horvath S. Tutorials for the WGCNA package \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html}
#'Fukushima A. (2013) DiffCorr: an R package to analyze and visualize differential correlations in biological networks. Gene, 10;518(1):209-14.
#'@export
#'@seealso \code{\link{WGCNA::cor}}, \code{\link{WGCNA::corPvalueStudent}}, \code{link{fetchDiffCorrNetwork}}, \code{\link{fetchGrinnNetwork}}, \pkg{\link{DiffCorr}}, \url{http://js.cytoscape.org/}
#'@examples
#'# Create metabolite-gene network from the list of metabolites using grinn ids and combine the grinn network to a differential correlation network of metabolites
#'kw <- c('X160','X300','X371')
#'dummyX1 <- rbind(nodetype=rep("metabolite"),mtcars[1:16,])
#'colnames(dummyX1) <- c('X1.1','X27967','X371','X4.1',letters[1:7])
#'rownames(dummyX1)[-1] <- paste0(rep("normal_"),1:16)
#'dummyX2 <- rbind(nodetype=rep("metabolite"),mtcars[17:32,])
#'colnames(dummyX2) <- c('X1.1','X27967','X371','X4.1',letters[1:7])
#'rownames(dummyX2)[-1] <- paste0(rep("cancer_"),1:16)
#'result <- fetchGrinnDiffCorrNetwork(txtInput=kw, from="metabolite", to="gene", datNormX1=dummyX1, datNormX2=dummyX2, pDiff=0.05)
#'plot(igraph::graph.data.frame(result$edges[,1:2], directed=F))
#'# Create metabolite-pathway network from the list of metabolites using grinn ids and combine the grinn network to a differential correlation network of metabolites and proteins
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
#'result <- fetchGrinnDiffCorrNetwork(txtInput=kw, from="metabolite", to="pathway", datNormX1=dummyX1, datNormX2=dummyX2, datNormY1=dummyY1, datNormY2=dummyY2, pDiff=0.05)
#'plot(igraph::graph.data.frame(result$edges[,1:2], directed=F))

fetchGrinnDiffCorrNetwork <- function(txtInput, from, to, filterSource=list(), returnAs="tab", dbXref="grinn", organism="'Homo sapiens'", 
                                  datNormX1=datNormX1, datNormX2=datNormX2, datNormY1=NULL, datNormY2=NULL, pDiff=1e-4, method="spearman"){
  basicnw = fetchGrinnNetwork(txtInput=txtInput,from=from,to=to,filterSource=filterSource,returnAs=returnAs,dbXref=dbXref,organism=organism)
  corrnw = fetchDiffCorrNetwork(datNormX1=datNormX1,datNormY1=datNormY1,datNormX2=datNormX2,datNormY2=datNormY2,pDiff=pDiff,method=method,returnAs=returnAs)
  if(nrow(basicnw$nodes)>0 && nrow(corrnw$nodes)>0){
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
  }else if(nrow(basicnw$nodes)>0 && nrow(corrnw$nodes)==0){
    cat("Formating and returning combined network ...\n")
    pair = basicnw$edges
    attb = basicnw$nodes
    cat("Found ",nrow(pair)," relationships...\n")
  }else if(nrow(basicnw$nodes)==0 && nrow(corrnw$nodes)>0){
    cat("Formating and returning combined network ...\n")
    pair = corrnw$edges
    attb = corrnw$nodes
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