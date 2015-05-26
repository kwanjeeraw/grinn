#'Compute an integrated network using information from grinn internal database
#'@description from the list of keywords, build an integrated network (grinn network) by connecting these keywords to a specified node type.
#'The keywords can be any of these node types: metabolite, protein, gene and pathway.
#'Grinn internal database contains the networks of the following types that can be quried: 
#'metabolite-protein, metabolite-protein-gene, metabolite-pathway, protein-gene, protein-pathway and gene-pathway.
#'@usage fetchGrinnNetwork(txtInput, from, to, filterSource, returnAs, dbXref, organism)
#'@param txtInput vector of keywords containing keyword ids e.g. txtInput = list('id1', 'id2'). 
#'The keyword ids are from the specified database, see \code{dbXref}. Default is grinn id e.g. X371.
#'@param from string of start node. It can be one of "metabolite","protein","gene","pathway".
#'@param to string of end node. It can be one of "metabolite","protein","gene","pathway".
#'@param filterSource string or list of pathway databases. The argument is required, if \code{from} or \code{to = "pathway"}, see \code{from} and \code{to}.
#'The argument value can be any of "SMPDB","KEGG","REACTOME" or combination of them e.g. list("KEGG","REACTOME").  
#'@param returnAs string of output type. Specify the type of the returned network. 
#'It can be one of "tab","json","cytoscape", default is "tab". "cytoscape" is the format used in Cytoscape.js
#'@param dbXref string of database name. Specify the database name used for the txtInput ids, see \code{txtInput}. 
#'It can be one of "grinn","chebi","kegg","pubchem","inchi","hmdb","smpdb","reactome","uniprot","ensembl","entrezgene". Default is "grinn".
#'If pubchem is used, it has to be pubchem SID (substance ID).
#'@param organism string of species in the following format: organism = "'species'". Default is "'Homo sapiens'".
#'@return list of nodes and edges. The list is with the following componens: edges and nodes. Return empty list if found nothing
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@export
#'@seealso \url{http://js.cytoscape.org/}
#'@examples
#'# Query genes by ENSEMBL ids and build a grinn network of gene-protein-metabolite
#'txtInput <- list('ENSG00000140459','ENSG00000143811','ENSG00000104524')
#'result <- fetchGrinnNetwork(txtInput, from="gene", to="metabolite", returnAs="tab", dbXref="ensembl")
#'library(igraph)
#'plot(graph.data.frame(result$edges[,1:2], directed=FALSE))
#'# Query metabolites by grinn ids and build a grinn network of metabolite-pathway
#'txtInput <- c('X371','X783')
#'result <- fetchGrinnNetwork(txtInput, from="metabolite", to="pathway", returnAs="json", organism="'Homo sapiens'")
#'# Query metabolites by grinn ids and build a network of metabolite-pathway using information from KEGG and REACTOME
#'txtInput <- c('X371','X783')
#'result <- fetchGrinnNetwork(txtInput, from="metabolite", to="pathway", filterSource=list("KEGG","REACTOME"), returnAs="tab", organism="'Homo sapiens'")
#'# Query proteins by uniprot ids and build a network of protein-pathway using information from SMPDB
#'txtInput <- list('P05108','Q53H96','P18463')
#'result <- fetchGrinnNetwork(txtInput, from="protein", to="pathway", filterSource="SMPDB", returnAs="cytoscape", dbXref="uniprot")

fetchGrinnNetwork <- function(txtInput, from, to, filterSource=list(), returnAs="tab", dbXref="grinn", organism="'Homo sapiens'"){ 
  out <- tryCatch(
  {
    tmparg <- try(from <- match.arg(from, c("metabolite","protein","gene","pathway"), several.ok = FALSE), silent = TRUE)
    if (class(tmparg) == "try-error") {
      stop("argument 'from' is not valid, choose one from the list: 'metabolite','protein','gene','pathway'")
    }
    tmparg <- try(to <- match.arg(to, c("pathway","protein","gene","metabolite"), several.ok = FALSE), silent = TRUE)
    if (class(tmparg) == "try-error") {
      stop("argument 'to' is not valid, choose one from the list: 'metabolite','protein','gene','pathway'")
    }
    tmparg <- try(dbXref <- match.arg(dbXref, c("grinn","chebi","kegg","pubchem","inchi","hmdb","smpdb","reactome","uniprot","ensembl","entrezgene"), several.ok = FALSE), silent = TRUE)
    if (class(tmparg) == "try-error") {
      stop("argument 'dbXref' is not valid, choose one from the list: 'grinn','chebi','kegg','pubchem','inchi','hmdb','smpdb','reactome','uniprot','ensembl','entrezgene'")
    }
    
    txtInput = unique(stringr::str_trim(as.character(txtInput))) #remove whiteline, duplicate
    
    #construct query string
    cat("Constructing query ...\n")
    ft = paste0(Hmisc::capitalize(from),Hmisc::capitalize(to))    
    len = length(txtInput)
    querystring = pathList[ft]
    if(tolower(dbXref) == 'grinn'){
      querystring = paste(querystring,'WHERE lower(from.GID) = lower(x)')
      txtInput = paste0("['",paste0(txtInput, collapse = "','"),"']")
    }else if(tolower(dbXref) == 'inchi'){
      querystring = paste(querystring,'WHERE from.InChI = x')
      txtInput = paste0("['",paste0(txtInput, collapse = "','"),"']")
    }else{
      querystring = paste(querystring,'WHERE ANY(y IN from.xref WHERE lower(y) = lower(x))')
      txtInput = paste0("['",paste0(dbXref,":", txtInput, collapse = "','"),"']") 
    }
    ##filter pathway sources
    if(from == "pathway" || to == "pathway"){filterSource=filterSource}else{filterSource=list()}
    if(length(filterSource)==0){
      querystring = paste(querystring,'RETURN DISTINCT ptw')
    }else{
      querystring = paste0(querystring,' AND (lower(rel.source) = lower(\'',paste0(filterSource, collapse = "\') OR lower(rel.source) = lower(\'"),'\')) RETURN DISTINCT ptw')
    }
    querystring = gsub("keyword", txtInput, querystring)
    querystring = gsub("species", organism, querystring)
  print(querystring)
    cat("Querying and returning network ...\n")
    network = curlRequestCypher(querystring)
    formatNetworkOutput(network,returnAs)
  },
  error=function(e) {
    message(e)
    cat("\n..No network found, RETURN empty list\n")
    cynetwork = list(nodes="", edges="")
    return(list(nodes=data.frame(),edges=data.frame())) # Choose a return value in case of error
  })    
  return(out)
}
#     #blockwise, split txtInput for shorter quried time
#     if(len>1000){#txtInput > 3900
#       cat("Long keywords, Partitioning quries ...\n")
#       bl = 500
#       tmpNode = data.frame()
#       x = floor(len/bl)
#       seqls = rep(1) + rep(seq(0,len-bl,bl), each=1)
#       for(i in 1:x){
#         sfr = seqls[i]
#         gto = (seqls[i]+bl)-1
#         if(i!=x){  
#           tmpkw = txtInput[sfr:gto]
#         } else{  
#           tmpkw = txtInput[sfr:len]
#         }
#         querystring = pathList[ft]
# 
#         cat("Converting and returning mapped ids ...\n")
#         tmpq = curlRequestCypher(querystring)
#         tmpNode = rbind(tmpNode,data.frame(gsub(paste0(dbXref,":"),"",tmpq)))
#         
#         if(tolower(dbXref) == 'grinn'){
#           querystring = paste(querystring,'WHERE lower(from.GID) = lower(x)')
#           tmpkw = paste0("['",paste0(tmpkw, collapse = "','"),"']")
#         }else if(tolower(dbXref) == 'inchi'){
#           querystring = paste(querystring,'WHERE from.InChI = x')
#           tmpkw = paste0("['",paste0(tmpkw, collapse = "','"),"']")
#         }else{
#           querystring = paste(querystring,'WHERE ANY(y IN from.xref WHERE lower(y) = lower(x))')
#           tmpkw = paste0("['",paste0(dbXref,":", tmpkw, collapse = "','"),"']") 
#         }
#         ##filter pathway sources
#         if(from == "pathway" || to == "pathway"){filterSource=filterSource}else{filterSource=list()}
#         if(length(filterSource)==0){
#           querystring = paste(querystring,'RETURN DISTINCT ptw')
#         }else{
#           querystring = paste0(querystring,' AND (lower(rel.source) = lower(\'',paste0(filterSource, collapse = "\') OR lower(rel.source) = lower(\'"),'\')) RETURN DISTINCT ptw')
#         }
#         querystring = gsub("keyword", tmpkw, querystring)
#         querystring = gsub("species", organism, querystring)
# print(querystring)
#         cat("Querying and returning network ...\n")
#         tmpq = curlRequestCypher(querystring)
#       }#end for
#       node = tmpNode
#       formatNetworkOutput(tmpq,returnAs)
#     }else{#txtInput less than 1000
#       querystring = pathList[ft]
#       if(tolower(dbXref) == 'grinn'){
#         querystring = paste(querystring,'WHERE lower(from.GID) = lower(x)')
#         txtInput = paste0("['",paste0(txtInput, collapse = "','"),"']")
#       }else if(tolower(dbXref) == 'inchi'){
#         querystring = paste(querystring,'WHERE from.InChI = x')
#         txtInput = paste0("['",paste0(txtInput, collapse = "','"),"']")
#       }else{
#         querystring = paste(querystring,'WHERE ANY(y IN from.xref WHERE lower(y) = lower(x))')
#         txtInput = paste0("['",paste0(dbXref,":", txtInput, collapse = "','"),"']") 
#       }
#       ##filter pathway sources
#       if(from == "pathway" || to == "pathway"){filterSource=filterSource}else{filterSource=list()}
#       if(length(filterSource)==0){
#         querystring = paste(querystring,'RETURN DISTINCT ptw')
#       }else{
#         querystring = paste0(querystring,' AND (lower(rel.source) = lower(\'',paste0(filterSource, collapse = "\') OR lower(rel.source) = lower(\'"),'\')) RETURN DISTINCT ptw')
#       }
#       querystring = gsub("keyword", txtInput, querystring)
#       querystring = gsub("species", organism, querystring)
# print(querystring)
#       cat("Querying and returning network ...\n")
#       network = curlRequestCypher(querystring)
#       formatNetworkOutput(network,returnAs)
#     }},