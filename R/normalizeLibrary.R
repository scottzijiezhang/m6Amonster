#' @title normalizeLibrary
#' @param readsOut Read counts data list from countReads() function.
#' @param X Study design or Grouping for the samples, should have 2 levels
#' @return norm_lib returns the list of normalized reads, study design and gene-bin names
#' @export
normalizeLibrary <- function(readsOut,X){

  input <- readsOut$reads[,1:length(readsOut$samplenames)]
  m6A <- readsOut$reads[,(1+length(readsOut$samplenames)):(2*length(readsOut$samplenames))]
  colnames(input) <- colnames(m6A) <- readsOut$samplenames

  ## check if geneBins already exist
  if( "geneBins" %in% names(readsOut) ){
    geneBins <- readsOut$geneBins
  }else{
    ## split gene and bin names
    aa <- strsplit(rownames(input), ",")
    gene.name <- unlist(lapply(aa, function(x){
      return(x[1])
    }))
    bin.name <- unlist(lapply(aa, function(x){
      return(x[2])
    }))
    geneBins <- data.frame(gene=gene.name,bin=as.integer(bin.name))
    rownames(geneBins) <- rownames(input)
  }

  ## Get input geneSum (gene level quantification)
  geneSum <- NULL
  for(i in 1:ncol(input) ){
    y <- input[,i]
    gene.sum <- tapply(y,gene.name,sum)
    geneSum <- cbind(geneSum,gene.sum)
  }
  colnames(geneSum) <- readsOut$samplenames

  size.input <- DESeq2::estimateSizeFactorsForMatrix(geneSum)

  norm.input <-t( t(input) / size.input )
  geneSum.norm <- t ( t(geneSum)/size.input)

  ## estimate enrichment using top IP count bins
  ave.ip <- rowMeans(m6A)
  ave.top <- order(ave.ip,decreasing = T)[1:round(0.01*length(ave.ip)[1])]

  geneCounts.window <- geneSum.norm[gene.name,]

  enrich <- as.data.frame(m6A[ave.top,]/geneCounts.window[ave.top,])
  enrich <- enrich[!apply(enrich,1, function(x){any(is.na(x)) | any(is.infinite(x))}),]

  size.enrich.deseq2 <- DESeq2::estimateSizeFactorsForMatrix(enrich[,1:length(X)])

  norm.ip <-t( t(m6A)/size.enrich.deseq2 )
  sizeFactor <- data.frame(input=size.input,ip=size.enrich.deseq2)
  ## check if geneBins already exist
  if( "geneBins" %in% names(readsOut) ){
    norm_lib <- c(readsOut, list('geneSum'=round(geneSum.norm),
                                 'norm.input'=norm.input,
                                 'norm.ip'=norm.ip,
                                 'sizeFactor'=sizeFactor,
                                 'X'=X)
    )
  }else{
    norm_lib <- c(readsOut, list('geneSum'=round(geneSum.norm),
                                 'norm.input'=norm.input,
                                 'norm.ip'=norm.ip,
                                 'sizeFactor'=sizeFactor,
                                 'geneBins'=geneBins,
                                 'X'=X)
    )
  }



  return(norm_lib)
}


#' @title normalizePeak
#' @param readsOut Read counts data list from countReads() function.
#' @param X Study design or Grouping for the samples, should have 2 levels
#' @export
normalizePeak <- function(readsOut,X){

  ## load data from input
  jointPeak_ip <- readsOut$jointPeak_ip
  jointPeak_input <- readsOut$jointPeak_input
  input <- readsOut$reads[,1:length(readsOut$samplenames)]
  colnames(input) <- readsOut$samplenames
  geneBins <- readsOut$geneBins

  ## Get input geneSum (gene level quantification)
  geneSum <- NULL
  for(i in 1:ncol(input) ){
    y <- input[,i]
    gene.sum <- tapply(y,geneBins$gene,sum)
    geneSum <- cbind(geneSum,gene.sum)
  }
  colnames(geneSum) <- readsOut$samplenames

  size.input <- DESeq2::estimateSizeFactorsForMatrix(geneSum)

  norm.jointPeak_input <-t( t(jointPeak_input) / size.input )
  geneSum.norm <- t ( t(geneSum)/size.input)

  ## Get the gene level input count for corresponding peaks
  geneCounts.peak <- geneSum.norm[geneBins[rownames(norm.jointPeak_input),"gene"],]
  enrich <- as.data.frame(jointPeak_ip/geneCounts.peak)
  enrich <- enrich[!apply(enrich,1, function(x){any(is.na(x)) | any(is.infinite(x))}),]

  size.enrich.deseq2 <- DESeq2::estimateSizeFactorsForMatrix(enrich[,1:length(X)])

  norm.jointPeak_ip <-t( t(jointPeak_ip)/size.enrich.deseq2 )
  sizeFactor <- data.frame(input=size.input,ip=size.enrich.deseq2)

  if("geneSum" %in% names(readsOut) ){
    data.out <- c(readsOut, list('norm.jointPeak_ip'=norm.jointPeak_ip))
  }else{
    data.out <- c(readsOut, list('geneSum'=round(geneSum.norm),
                               'norm.jointPeak_ip'=norm.jointPeak_ip,
                               'sizeFactor'=sizeFactor,
                               'X'=X))
  }

  return(data.out)
}
