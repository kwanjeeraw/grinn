#' \code{fetchNode} query for an entity or list of entities from grinn internal database
#'@description from the list of keywords, query for entities from the grinn internal database. Quried results include entity information and relationship information.
#'The keywords can be any of these node types: metabolite, protein, gene and pathway.
#'Relationship information includes information of adjacent nodes, name and additional attributes of relationships.
#'@usage fetchNode(txtInput, nodetype, searchField, exactMatch, returnAs, dbXref, organism)
#'@param txtInput list of keywords, can be one of grinn id, name, synonym, database xref, inchi e.g. txtInput = list(name1, name2). 
#'The type of keywords needs to be provided in the argument \code{searchField}, see \code{searchField}. Default is grinn id e.g. X371.
#'@param nodetype string of entity type. It can be one of "metabolite","protein","gene","pathway".
#'@param searchField string of keyword property to search for, see \code{txtInput}. It can be one of "grinn","name","synonym","xref","inchi". Default is "grinn".
#'If \code{searchField} = "xref", the argument \code{dbXref} is required, see \code{dbXref}.
#'If \code{searchField} = "name" or "synonym", the argument \code{exactMatch} is required, see \code{exactMatch}.
#'@param exactMatch boolean value. The argument is required if \code{searchField} = "name" or "synonym", see \code{searchField}. Default is TRUE.
#'@param returnAs string of output type. It can be one of "list","json", default is "list".
#'@param dbXref string of database name. Specify the database name used for the txtInput when the argument \code{searchField}="xref", see \code{txtInput} and \code{searchField}. 
#'It can be one of "grinn",chebi","kegg","pubchem","hmdb","smpdb","reactome","uniprot","ensembl","entrezgene". Default is "grinn".
#'If pubchem is used, it has to be pubchem SID (substance ID).
#'@param organism string of species in the following format: organism = "'species'". Default is "'Homo sapiens'". \code{organism} required, if \code{nodetype}="metabolite".
#'@return list of nodes and first neighborhoods. Return empty list if found nothing.
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@export
#'@examples
#'# Query metabolites by grinn ids
#'txtInput <- list('X371','X783')
#'result <- fetchNode(txtInput, nodetype="metabolite")
#'# Query genes by KEGG ids
#'txtInput <- list('hsa:4514','hsa:4537')
#'result <- fetchNode(txtInput, nodetype="gene", searchField="xref", dbXref="kegg", returnAs="list", organism="'Homo sapiens'")
#'# Query proteins by names
#'txtInput <- list('14-3-3 protein beta/alpha','ARL14 effector protein-like','6-phosphogluconolactonase')
#'result <- fetchNode(txtInput, nodetype="protein", searchField="name", returnAs="list", organism="'Homo sapiens'")

fetchNode <- function(txtInput, nodetype, searchField="grinn", exactMatch=TRUE, returnAs="list", dbXref="grinn", organism="'Homo sapiens'"){ 
  out <- tryCatch(
  {

    tmparg <- try(nodetype <- match.arg(nodetype, c("metabolite","protein","gene","pathway"), several.ok = FALSE), silent = TRUE)
    if (class(tmparg) == "try-error") {
      stop("argument 'nodetype' is not valid, choose one from the list: 'metabolite','protein','gene','pathway'")
    }
    tmparg <- try(searchField <- match.arg(searchField, c("grinn","name","synonym","xref","inchi"), several.ok = FALSE), silent = TRUE)
    if (class(tmparg) == "try-error") {
      stop("argument 'searchField' is not valid, choose one from the list: 'grinn','name','synonym','xref','inchi'")
    }
    tmparg <- try(dbXref <- match.arg(dbXref, c("grinn","chebi","kegg","pubchem","inchi","hmdb","smpdb","reactome","uniprot","ensembl","entrezgene"), several.ok = FALSE), silent = TRUE)
    if (class(tmparg) == "try-error") {
      stop("argument 'dbXref' is not valid, choose one from the list: 'grinn','chebi','kegg','pubchem','inchi','hmdb','smpdb','reactome','uniprot','ensembl','entrezgene'")
    }
    
    txtInput = unique(stringr::str_trim(as.character(txtInput))) #remove whiteline, duplicate
    #construct query string
    nodetype=Hmisc::capitalize(nodetype)
    if(searchField == 'grinn'||searchField == 'xref'||searchField == 'inchi'){exactMatch=TRUE}
    
    isString =  searchField %in% propertyList$stringVal
    if (exactMatch == TRUE && isString == TRUE) {
      querystring = nodeList["exactMatch"]
    }else if (exactMatch == FALSE && isString == TRUE) {
      querystring = nodeList["regexMatch"] 
    }else if (exactMatch == TRUE && isString == FALSE) {
      querystring = nodeList["exactCollection"] 
    }else{
      querystring = nodeList["regexCollection"]
    }

    if(searchField == 'grinn'){
      querystring = gsub("property", "GID", querystring)
      txtInput = paste0("['",paste0(txtInput, collapse = "','"),"']")
    }else if(searchField == 'xref'){
      querystring = gsub("property", searchField, querystring)
      txtInput = paste0("['",paste0(dbXref,":", txtInput, collapse = "','"),"']") 
    }else{
      querystring = gsub("property", searchField, querystring)
      txtInput = paste0("['",paste0(txtInput, collapse = "','"),"']")
    }
    
    querystring = gsub("keyword", txtInput, querystring)
    querystring = gsub("label", nodetype, querystring)
    if(nodetype!="Metabolite"){
      querystring = paste(querystring,'AND node.organism =',organism)
    }
    querystring = paste(querystring,'RETURN DISTINCT node')  
print(querystring)     
    node = curlRequestCypher(querystring)
    formatNodeOutput(node,returnAs)
  },
  error = function(e) {
    message(e)
    cat("\n..No match, RETURN empty list")
    return(list()) # Choose a return value in case of error
  })    
  return(out)
}