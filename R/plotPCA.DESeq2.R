plotPCA.DESeq2 <- function(data,group = NULL, returnPC = FALSE){
  if(is.null(group)){
    colData <- data.frame(group=colnames(data))
  }else{
    colData <- data.frame(group=group)
  }
  
  rownames(colData) <- colnames(data)
  dds <- DESeq2::DESeqDataSetFromMatrix(data,colData,design = ~group)
  cat("Using regularized log transformation of DESeq2 to tranform data...\n")
  rld <- DESeq2::rlog(dds)
  cat("Plot PCA using the rlog transformed data...\n")
  DESeq2::plotPCA(rld,intgroup = "group")
  if(returnPC){
    PCs <- DESeq2::plotPCA(rld,intgroup = "group",returnData=TRUE)
  }
  return(PCs)
}