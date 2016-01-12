#'Combine two networks
#'@description The function combines two networks at a time. The networks might be computed from grinn functions such as \code{fetchGrinnNetwork}, \code{fetchCorrGrinnNetwork}, 
#'\code{fetchDiffCorrGrinnNetwork}, or from users' provided  networks. 
#'The provided network must be in the following format:
#'\code{nodes} = data.frame(id,nodename,nodetype) #nodes contain at least 3 columns
#'\code{edges} = data.frame(source,target,relname) #edges contain at least 3 columns
#'\code{nwX} or \code{nwY} = list(nodes,edges)
#'@usage combineNetwork(nwX, nwY, returnAs)
#'@param nwX list of nodes and edges to be combined.
#'@param nwY list of nodes and edges to be combined.
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape", default is "tab". "cytoscape" is the format used in Cytoscape.js
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Return empty list if error.
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@export
#'@seealso \url{http://js.cytoscape.org/}, \code{\link{fetchGrinnNetwork}}, \code{\link{fetchCorrNetwork}}, \code{\link{fetchWGCNAModule}}, \code{\link{fetchDiffCorrNetwork}}
#'@examples
#'# Create metabolite-protein network from the list of metabolites using grinn ids and combine the grinn network to a correlation network
#'kw <- c('G160','G300','G371')
#'grinnnw <- fetchGrinnNetwork(txtInput=kw, from="metabolite", to="protein")
#'# Compute a correlation network of metabolites and expand to a grinn network of metabolite-protein
#'dummy <- rbind(nodetype=rep("metabolite"),t(mtcars))
#'colnames(dummy) <- c('G1.1','G15603','G371','G17295',paste0('G',sample(400:22000, 28)))
#'corrnw <- fetchCorrGrinnNetwork(datX=dummy, corrCoef=0.7, pval=1e-12, method="spearman", returnAs="tab", xTo="protein")
#'result <- combineNetwork(grinnnw,corrnw)
#'library(igraph)
#'plot(graph.data.frame(result$edges[,1:2], directed=FALSE))
#'# Create metabolite-protein network from the list of metabolites using grinn ids and combine the grinn network to a provided network
#'txtInput <- list('G371','G783','G1.1')
#'grinnnw <- fetchGrinnNetwork(txtInput, from="metabolite", to="protein")
#'nwY = list()
#'nwY$nodes = data.frame(id=c("G371","G783","XXX","YYY"),nodename=c("G371","G783","XXX","YYY"),nodetype=c("metabolite","metabolite","protein","protein"),stringsAsFactors = FALSE)
#'nwY$edges = data.frame(source=c("G371","G783"),target=c("XXX","YYY"),corr_coef=c(-0.368,0.385),pval=c(0.000927,0.000497),reltype=c("Metabolite_Protein","Metabolite_Protein"),relname=c("CORRELATION","CORRELATION"),stringsAsFactors = FALSE)
#'result <- combineNetwork(grinnnw,nwY)

combineNetwork <- function(nwX, nwY, returnAs="tab"){
  tryCatch({
    if(!(length(colnames(nwY$nodes)) >= 3) || !(length(colnames(nwY$edges)) >= 4)){
      stop("incorrect number of columns or incorrect column names, \n
           require at least 'id, nodename, nodetype' for nodelist and 'source, target, reltype, relname' for edgelist")
    }
    
    cat("Formating and returning combined network ...\n")
    
    attb = plyr::rbind.fill(nwX$nodes,nwY$nodes)
    pair = plyr::rbind.fill(nwX$edges,nwY$edges)
    pair = pair[!duplicated(pair[c("source","target")]), ] #remove duplicate edges
    attb = attb[!duplicated(attb[,1]),] #remove duplicate nodes
    colbs = colnames(nwX$nodes)
    colot = colnames(nwY$nodes)
    cat("Found ",nrow(pair)," relationships...\n")
    
    out = switch(returnAs,
                 tab = list(nodes=attb, edges=pair),
                 json = list(nodes=jsonlite::toJSON(attb), edges=jsonlite::toJSON(pair)),
                 cytoscape = createCyNetwork(attb, pair),
                 stop("incorrect return type"))
  }, error = function(err) {
    print(err)
    pair = data.frame()
    attb = data.frame()
    out <- list(nodes=attb, edges=pair) #return empty list
  }) # END tryCatch
}