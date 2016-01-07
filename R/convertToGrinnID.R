#'Convert other database IDs to grinn IDs
#'@description Convert the list of ids to grinn IDs
#'grinn IDs are recommended to be used in several functions including \code{combineNetwork}, \code{fetchGrinnCorrNetwork}, \code{fetchGrinnModuNetwork},
#'\code{fetchGrinnDiffCorrNetwork}, \code{fetchCorrGrinnNetwork}, \code{fetchModuGrinnNetwork}, \code{fetchDiffCorrGrinnNetwork}. 
#'@usage convertToGrinnID(txtInput, nodetype, dbXref)
#'@param txtInput list of keywords containing ids e.g. \code{txtInput} = list('id1', 'id2'). 
#'The ids are from a database, see \code{dbXref}.
#'@param nodetype string of entity type. It can be one of "metabolite","protein","gene","pathway".
#'@param dbXref string of database name. Specify the database name used for the txtInput, 
#'see \code{txtInput}. It can be one of "chebi","kegg","pubchem","inchi","hmdb","smpdb","reactome","uniprot","ensembl","entrezgene".
#'@return data frame of ids. Return empty list if found nothing.
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@export
#'@examples
#'# Convert genes from ENSEMBL ids to grinn ids
#'txtInput <- list('ENSG00000140459','ENSG00000143811','ENSG00000104524')
#'result <- convertToGrinnID(txtInput, nodetype="gene", dbXref="ensembl")
#'# Convert metabolites from chebi ids to Grinn ids
#'txtInput <- c(371,783)
#'result <- convertToGrinnID(txtInput, nodetype="metabolite", dbXref="chebi")
#'# Convert metabolites from kegg ids to Grinn ids, read input from csv or txt file
#'library(grinn)
#'data(keggids)
#'result <- convertToGrinnID(unlist(keggids), nodetype="metabolite", dbXref="kegg")

convertToGrinnID <- function(txtInput, nodetype, dbXref){
  tmparg <- try(nodetype <- match.arg(nodetype, c("metabolite","protein","gene","pathway"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'nodetype' is not valid, choose one from the list: 'metabolite','protein','gene','pathway'")
  }
  tmparg <- try(dbXref <- match.arg(dbXref, c("chebi","kegg","pubchem","inchi","hmdb","smpdb","reactome","uniprot","ensembl","entrezgene"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'dbXref' is not valid, choose one from the list: 'chebi','kegg','pubchem','inchi','hmdb','smpdb','reactome','uniprot','ensembl','entrezgene'")
  }
  #construct query string
  nodetype=Hmisc::capitalize(nodetype)
  txtInput = unique(stringr::str_trim(as.character(txtInput))) #remove whiteline, duplicate
  len = length(txtInput)
  #blockwise, split txtInput for shorter quried time
  if(len>1000){#txtInput > 1000
    cat("Long keywords, Partitioning quries ...\n")
    bl = 500
    tmpNode = data.frame()
    x = floor(len/bl)
    seqls = rep(1) + rep(seq(0,len-bl,bl), each=1)
    for(i in 1:x){
      sfr = seqls[i]
      gto = (seqls[i]+bl)-1
      if(i!=x){  
        tmpkw = txtInput[sfr:gto]
      } else{  
        tmpkw = txtInput[sfr:len]
      }
      if (tolower(dbXref) == 'inchi') {
        querystring = nodeList["exactMatch"]
        querystring = gsub("property", "inchi", querystring)
        tmpkw = paste0("['",paste0(tmpkw, collapse = "','"),"']")
      }else{
        querystring = nodeList["exactCollection"]
        querystring = gsub("property", "xref", querystring)
        tmpkw = paste0("['",paste0(dbXref,":", tmpkw, collapse = "','"),"']") 
      }
      querystring = gsub("label", nodetype, querystring)
      querystring = paste(querystring,'RETURN DISTINCT x, node.GID')
      querystring = gsub("keyword", tmpkw, querystring)
print(querystring)
      cat("Converting and returning mapped ids ...\n")
      tmpq = curlRequestCypher(querystring)
      tmpNode = rbind(tmpNode,data.frame(gsub(paste0(dbXref,":"),"",tmpq)))
    }
    node = tmpNode
    colnames(node) = c(paste0("FROM_",dbXref),"GRINNID")
    out <- node
  }else{#txtInput less than 1000
    if (tolower(dbXref) == 'inchi') {
      querystring = nodeList["exactMatch"]
      querystring = gsub("property", "inchi", querystring)
      txtInput = paste0("['",paste0(txtInput, collapse = "','"),"']")
    }else{
      querystring = nodeList["exactCollection"]
      querystring = gsub("property", "xref", querystring)
      txtInput = paste0("['",paste0(dbXref,":", txtInput, collapse = "','"),"']") 
    }
    querystring = gsub("label", nodetype, querystring)
    querystring = paste(querystring,'RETURN DISTINCT x, node.GID')
    querystring = gsub("keyword", txtInput, querystring)
print(querystring)
    cat("Converting and returning mapped ids ...\n")
    node = curlRequestCypher(querystring)
    colnames(node) = c(paste0("FROM_",dbXref),"GRINNID")
    out <- data.frame(gsub(paste0(dbXref,":"),"",node))
  }
}