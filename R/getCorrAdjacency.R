#'compute correlation coefficients, pvalues and directions
#'@description compute correlation coefficients, pvalues and relation directions from normalized data 
#'using \code{WGCNA::cor} and \code{WGCNA::corPvalueStudent} to return outputs for further uses by 
#'\code{fetchCorrNetwork}, \code{fetchDiffCorrNetwork}.
#'datNormX and datNormY are matrices in which rows are samples and columns are entities, require the first row to specify entity types e.g. metabolite.
#'If datNormY is given, then the correlations between the columns of datNormX and the columns of datNormY are computed.
#'Otherwise the correlations of the columns of datNormX are computed.
#'@references Langfelder P. and Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9:559 
#'@references Dudoit S., Yang YH., Callow MJ. and Speed TP. (2002) Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA, 12:111
#'@references Langfelder P. and Horvath S. Tutorials for the WGCNA package \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html}
#'@seealso \code{\link{cor}}, \code{\link{corPvalueStudent}}, \code{\link{fetchCorrNetwork}}, \code{\link{fetchDiffCorrNetwork}}
#'datMet = read.delim("~/grinn_sample/Lung_MET.txt", header=TRUE, stringsAsFactors=FALSE)
#'datPrt = read.delim("~/grinn_sample/Lung_PRT.txt", header=TRUE, stringsAsFactors=FALSE)
#'Convert kegg ids to grinn ids
#'grinnID = convertToGrinnID(txtInput=datMet$KEGG.id, nodetype="metabolite", dbXref="kegg")
#'grinnID = grinnID[!duplicated(grinnID[,1]),] #keep the first mapped id
#'internal function, called by apply for formatting ids
#'mapToInput = function(x){
#'  id = which(grinnID$FROM_kegg == x[2])
#'  out = ifelse(length(id)>0,as.character(grinnID$GRINNID[id]),x[1]) }
#'formatting input data
#'datMet = t(datMet)
#'colnames(datMet) = apply(datMet,2,mapToInput)
#'datMet = datMet[-c(1,2),]
#'formatting input data
#'datPrt = t(datPrt)
#'colnames(datPrt) = datPrt[1,]
#'datPrt = datPrt[-1,]
#'corAdj = getCorrAdjacency(datNormX=datMet, datNormY=datPrt, method="spearman")

getCorrAdjacency <- function (datNormX, datNormY, method) 
{
  if(!("nodetype" %in% rownames(datNormX))){
    stop("can't define type of entities, missing nodetype")
  }
  tmparg <- try(method <- match.arg(method, c("pearson","kendall","spearman"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'method' is not valid, choose one from the list: 'pearson', 'kendall', 'spearman'")
  }
  colnames(datNormX) = paste0(colnames(datNormX),"_TYPE_",datNormX[1,])
  datNormX = apply(datNormX[-1,],2,as.numeric)
  if (!is.null(datNormY)) {
    if(!("nodetype" %in% rownames(datNormY))){
      stop("can't define type of entities, missing nodetype")
    }
    colnames(datNormY) = paste0(colnames(datNormY),"_TYPE_",datNormY[1,])
    datNormY = apply(datNormY[-1,],2,as.numeric)

    cat("Computing correlation coefficients ...\n")
    cor_mat = WGCNA::cor(datNormX,datNormY,method=method)
    corP_mat = WGCNA::corPvalueStudent(cor_mat,nrow(datNormX))    
    #format output
    cor_df = reshape2::melt(cor_mat) #matrix to data frame
    colnames(cor_df) = c("source","target","corr_coef")
    corP_df = reshape2::melt(corP_mat) #matrix to data frame
    colnames(corP_df) = c("source","target","pval")
    #corAdj_df = p.adjust(p=corP_df$pval,method="none") #p-adjusted, no adjustment in this version
    corSign_df = sign(cor_df$corr_coef) #keep direction, -1 negative, 1 positive
    edgeData = cbind(cor_df,pval=corP_df$pval,direction=corSign_df)
  }
  if (is.null(datNormY)) {#Assign values to a square matrix
    cat("Computing correlation coefficients ...\n")
    cor_mat = WGCNA::cor(datNormX,method=method)
    corP_mat = WGCNA::corPvalueStudent(cor_mat,nrow(datNormX))
    nRow = nrow(cor_mat)
    nNames = dimnames(cor_mat)[[1]]
    rowMat = matrix(c(1:nRow), nRow, nRow, byrow = TRUE)
    colMat = matrix(c(1:nRow), nRow, nRow)
    dstRows = as.dist(rowMat)
    dstCols = as.dist(colMat)
    edgeData = data.frame(fromNode = nNames[dstRows], toNode = nNames[dstCols], corr_coef = cor_mat[lower.tri(cor_mat)], 
                            pval = corP_mat[lower.tri(corP_mat)], direction = sign(cor_mat[lower.tri(cor_mat)]))  
    colnames(edgeData) = c("source","target","corr_coef", "pval", "direction")
  }
  #out <- cbind(cor_df,pval=corP_df$pval,adjPval=corAdj_df,direction=corSign_df) #p-adjusted, no adjustment in this version
  out <- edgeData
}