#' \code{formatWGCNAModule} format resulting WGCNA module
#'@description Modify the R codes from WGCNA tutorials, to return module output for further uses by 
#'\code{fetchGrinnModuNetwork}, \code{fetchModuGrinnNetwork}, \code{fetchWGCNAModule}.
#'@references 
#'Langfelder P. and Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9:559 
#'Dudoit S., Yang YH., Callow MJ. and Speed TP. (2002) Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA, 12:111
#'Langfelder P. and Horvath S. Tutorials for the WGCNA package \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html}
#'@seealso \code{\link{exportNetworkToCytoscape}}, \code{\link{TOMsimilarityFromExpr}}, \code{fetchGrinnModuNetwork}, \code{fetchModuGrinnNetwork}, \code{fetchWGCNAModule}

formatWGCNAModule <- function(datNorm, softPower, mdColors, threshold, nodetype, returnAs)
{ 
  #Let users enter module color(s) they want at the terminal.
  cat("Waiting for user response ...\n")
  md <- readline(prompt="Enter module color(s), please separate each color by space: ")
  modules = unlist(strsplit(md, '\\s+'))
  cat("Generating output from the selected module(s) ...\n")
  #Select entities of wanted module(s)
  entities = colnames(datNorm)
  inModule = is.finite(match(mdColors, modules))
  modEntities = entities[inModule]
  #Calculate topological overlap
  TOM = WGCNA::TOMsimilarityFromExpr(datNorm, power=softPower)
  #Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modEntities, modEntities)
  #Call WGCNA function to get edge and node list
  options(stringsAsFactors = FALSE)
  cyt = WGCNA::exportNetworkToCytoscape(modTOM, threshold = threshold,
                                 nodeNames = modEntities,nodeAttr = mdColors[inModule])
  if(nrow(cyt$nodeData) > 0){
    cat("Formating and returning network module ...\n")
    pair = cyt$edgeData[,1:3]
    pair$relname = "MODULE"
    colnames(pair) = c("source","target","corr_coef","relname")
    attb = cyt$nodeData[,c(1,3)]
    attb$id = attb$nodeName
    attb$nodetype = Hmisc::capitalize(nodetype)
    attb = attb[,c(3,1,2,4)]
    colnames(attb) = c("id","nodename","modulecolor","nodetype")
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