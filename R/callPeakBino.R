#' @title callPeakFisher
#' @param readsOut The list object of countReads function.
#' @param min_counts The minimal number of reads present in a bin to be called a peak.
#' @param peak_cutoff_fdr The cutoff of fdr of fisher's exact test to call peak.
#' @param peak_cutoff_oddRatio The minimal oddRatio of fisher's exact test to call peak.
#' @param threads The number of threads to use. Default uses 1 threads.
#' @export
callPeakBinomial <- function(readsOut, min_counts = 15, peak_cutoff_fdr = 0.05 , peak_cutoff_oddRatio = 1, threads = 1){
  input <- as.matrix(readsOut$reads[,1:length(readsOut$samplenames)])
  m6A <- as.matrix(readsOut$reads[,(1+length(readsOut$samplenames)):(2*length(readsOut$samplenames))])
  T0 <- colSums(input)
  T1 <- colSums(m6A)
  colnames(input) <- colnames(m6A) <- readsOut$samplenames

  registerDoParallel( min( length( readsOut$samplenames) ,threads) )

  bino_pvalue <- foreach(j = 1:length( readsOut$samplenames), .combine = cbind )%dopar%{
     .Bino_test( m6A[,j], input[,j], T1[j], T0[j])
  }


  oddRatio <-  t(t(m6A)/T1) / t(t(input)/T0)

  ## total read count threshold
  above_thresh_counts <- ( (input + m6A ) >= min_counts )

  fdr <- foreach(j = 1:length( readsOut$samplenames), .combine = cbind )%dopar%{
    tmpFDR <- p.adjust( bino_pvalue[,j] , method = "fdr")
    ID <- which(above_thresh_counts[,j])
    tmpFDR[ ID ] <-  p.adjust( bino_pvalue[ID,j] , method = "fdr")
    tmpFDR
  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)

  test_peak <- (fdr < peak_cutoff_fdr & oddRatio > peak_cutoff_oddRatio &  above_thresh_counts )

  colnames(test_peak) <- readsOut$samplenames
  rownames(test_peak) <- rownames(input)

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

  ## check if geneBins already exist
  if( "geneBins" %in% names(readsOut) ){
    data.out <- c(readsOut,list('peakCallResult' = test_peak) )
  }else{
    data.out <- c(readsOut,list('geneBins'=geneBins,
                                'peakCallResult' = test_peak) )
  }

  return(data.out)

}


#' @title ctest
#' @return data frame with p-values and odds ratio
.Bino_test <- function(IP, input, IP_overall, input_overall, pseudo_count=1){

  IP <- pmax(IP,pseudo_count)
  input <- pmax(input,pseudo_count)

  p <- IP_overall/(IP_overall+input_overall)

  pvalue <-  pbinom(IP-1, size = (IP+input), p, lower.tail = FALSE )

  return(pvalue)
}
