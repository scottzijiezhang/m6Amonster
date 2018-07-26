#' @title getTPM
#' @param x The m6Amonster object
#' @param meanFragmentLength The mean length of RNA fragment (insert of RNA library). Default is 150bp.
#' @param  normalize Logical indicating whether normalized TPM or raw TPM should be returned.
#' @export
getTPM <- function(x, meanFragmentLength = 150,normalize = T){

  input <- x$reads[,1:length(x$samplenames)]
  gene.name <- x$geneBins$gene
  geneSum <- apply(input,2,function(y){
    tapply(y,gene.name,sum)
  })
  colnames(geneSum) <- x$samplenames

  genes <- rownames(x$geneSum)

  cat("calculating gene length...\n")
  geneLength <- sapply(genes,function(yy){
    sum( as.data.frame( x$geneModel[[yy]] )$width )
  })


  effLength <- geneLength - meanFragmentLength
  effLength <- sapply(effLength,max,1) ## remove effective length smaller than 0.

  cat("computing TPM from read counts...\n")

  rate <- apply(geneSum,2,function(aa){aa/effLength})
  totalCounts <-colSums(rate)

  tpm <- t( t(rate)/totalCounts ) *1e6

  size.tpm <- DESeq2::estimateSizeFactorsForMatrix(tpm)
  tpm_norm <- t(t(tpm)/ size.tpm)

  if(normalize){
    return(tpm_norm)
  }else{
    return(tpm)
  }
}
