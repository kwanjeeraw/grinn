#'Compute a network module correlated to a phenotypic feature
#'@description identify correlation between the input omics data e.g. normalized gene expression data, and phenotypic data e.g. weight.
#'The function wraps around important aspects of \pkg{WGCNA} including \code{blockwiseModules}, \code{cor}, \code{corPvalueStudent}, \code{labeledHeatmap}.
#'These aspects automatically perform correlation network construction, module detection, and display module-phenotype correlations.
#'A module or the combination of modules can be selected from the heatmap of module-phenotype correlations for including in the network output, see more details below. 
#'@usage fetchWGCNAModule(datX, datPheno, sfPower, minModuleSize, threshold, returnAs)
#'@param datX data frame containing normalized, quantified omics data e.g. expression data, metabolite intensities. 
#'Columns correspond to entities e.g. genes, metabolites, and rows to samples e.g. normals, tumors. Require 'nodetype' at the first row to indicate the type of entities in each column.
#'@param datPheno data frame containing phenotypic data e.g. weight, age, insulin sensitivity. Columns correspond to phenotypes and rows to samples e.g. normals, tumors. 
#'@param sfPower numerical value of soft-thresholding power for correlation network construction. It is automatically estimated using \code{pickSoftThreshold}, or it can be defined by users.
#'@param minModuleSize numerical value of minimum module size for module detection.
#'@param threshold numerical value to define the minimum value of similarity threshold, from 0 to 1, to include edges in the output.
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape". "cytoscape" is the format used in Cytoscape.js
#'@details
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
#'@seealso \code{\link{cor}}, \code{\link{corPvalueStudent}}, \code{\link{pickSoftThreshold}}, \code{\link{blockwiseModules}}, \code{\link{labeledHeatmap}}, \code{\link{exportNetworkToCytoscape}}, \url{http://js.cytoscape.org/}
#'@examples
#'# Compute a correlation of metabolite module to phenotypic data
#'#library(grinn)
#'#data(dummy)
#'#data(dummyPheno)
#'#result <- fetchWGCNAModule(datX=dummy, datPheno=dummyPheno, minModuleSize=5, threshold=0.2)
#'# enter module color(s) seperate by space:yellow brown purple
#'#library(igraph)
#'#plot(graph.data.frame(result$edges[,1:2], directed=FALSE))

fetchWGCNAModule <- function(datX, datPheno, sfPower=NULL, minModuleSize = 10, threshold = 0.5, returnAs="tab"){
  if(!("nodetype" %in% rownames(datX))){
    stop("can't define type of entities, missing nodetype")
  }
  nodetype = datX[1,1] #get nodetype
  datX = apply(datX[-1,],2,as.numeric) #remove unwant row and transform to numerical matrix
  #Call WGCNA function to choose power
  if(is.null(sfPower)){
    cat("Calculating and returning soft-thresholding power ...\n")
    #Define a set of soft-threshold powers
    powers = 1:20
    #Use the network topology analysis function
    sft = WGCNA::pickSoftThreshold(datX, powerVector = powers, verbose = 5)
    #Choose the lowest power where the scale-free topology fit index ~ 0.90
    ind = which((-sign(sft$fitIndices[,3])*sft$fitIndices[,2]) >= 0.9)
    softPower = sft$fitIndices[,1][ind[1]]
  }else{#Otherwise use users defined power 
    softPower = sfPower
  }
  #Automatic correlation network construction and module detection with WGCNA function
  cat("Building a network and detecting network modules ...\n")
  corrMD = WGCNA::blockwiseModules(datX, power = softPower, corType = "pearson", minModuleSize = minModuleSize, verbose = 3)
  mdLabels = corrMD$colors
  mdColors = WGCNA::labels2colors(corrMD$colors)
  MEs = corrMD$MEs
  nSamples = nrow(datX)
  mdTraitCor = WGCNA::cor(MEs, datPheno, use = "p")
  mdTraitPvalue = WGCNA::corPvalueStudent(mdTraitCor, nSamples)
  textMatrix = paste(signif(mdTraitCor, 2), "\n(",signif(mdTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(mdTraitCor)
  #Show a heatmap of modules-features correlations and pvalues with WGCNA function
  cat("Displaying a heatmap of modules-phenotypes correlations ...\n")
  par(mar = c(6, 8.5, 3, 3))
  WGCNA::labeledHeatmap(Matrix = mdTraitCor, xLabels = colnames(datPheno), yLabels = colnames(MEs), ySymbols = names(MEs), colorLabels = FALSE,
                 colors = WGCNA::blueWhiteRed(60), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.6, zlim = c(-1,1), main = paste("Module-phenotype correlations"))
  out <- formatWGCNAModule(datX=datX,softPower=softPower,mdColors=mdColors,threshold=threshold,nodetype=nodetype,returnAs=returnAs)  
}