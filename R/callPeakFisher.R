#' @title callPeakFisher
#' @param readsOut The list object of countReads function.
#' @param min_counts The minimal number of reads present in a bin to be called a peak.
#' @param peak_cutoff_fdr The cutoff of fdr of fisher's exact test to call peak.
#' @param peak_cutoff_oddRatio The minimal oddRatio of fisher's exact test to call peak.
#' @param threads The number of threads to use. Default uses 1 threads.
#' @export
callPeakFisher <- function(readsOut, min_counts = 15, peak_cutoff_fdr = 0.05 , peak_cutoff_oddRatio = 1, threads = 1){
  input <- as.matrix(readsOut$reads[,1:length(readsOut$samplenames)])
  m6A <- as.matrix(readsOut$reads[,(1+length(readsOut$samplenames)):(2*length(readsOut$samplenames))])
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

  ## number of genes to call peaks
  batch_id_list <- unique(geneBins$gene)
  num_batch_ids <- length(batch_id_list)
  cat("Calling peaks for ",num_batch_ids, "genes... \n")
  registerDoParallel( cores = threads)
  start_time <- Sys.time()
  cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
  cat(paste("Using",getDoParWorkers(),"thread(s) to call peaks in continuous bins...\n"))
  peak_call_batches <- foreach( i = 1:num_batch_ids, .combine = rbind)%dopar%{

    idx_batch <- which(geneBins$gene == batch_id_list[i])
    batch_input <- input[idx_batch,]
    batch_m6A <- m6A[idx_batch,]
    overall_input <- round( apply(batch_input, 2, median, na.rm = TRUE)  )
    overall_m6A <- round( apply(batch_m6A, 2, median, na.rm = TRUE)  )

    ## loop through all sample
    fisher_exact_test_p <- NULL
    fisher_exact_test_oddRatio <- NULL
    for(j in 1:length(overall_input) ){
      fisher_result <- t( mapply(.fisher_exact_test, batch_m6A[,j], batch_input[,j], overall_m6A[j], overall_input[j]) )
      fisher_exact_test_p <- cbind(fisher_exact_test_p,fisher_result[,1])
      fisher_exact_test_oddRatio <- cbind(fisher_exact_test_oddRatio,fisher_result[,2] )
    }

    above_thresh_counts <- ( (batch_input + batch_m6A) >= min_counts )

    fisher_exact_test_fdr <- matrix(1,nrow = nrow(fisher_exact_test_p),ncol = ncol(fisher_exact_test_p))
    if(sum(rowSums(above_thresh_counts)> (length(overall_input)/2))>1){
      fisher_exact_test_fdr[rowSums(above_thresh_counts)> (length(overall_input)/2) ,] <- apply(fisher_exact_test_p[which(rowSums(above_thresh_counts)> (length(overall_input)/2)) ,] , 2, p.adjust, method = 'fdr')
    }

    fisher_exact_test_peak <- (fisher_exact_test_fdr < peak_cutoff_fdr &
                                 fisher_exact_test_oddRatio > peak_cutoff_oddRatio &
                                 above_thresh_counts)

    fisher_exact_test_peak
  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
  end_time <- Sys.time()
  cat(paste("Time used to call peaks:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
  colnames(peak_call_batches) <- readsOut$samplenames

  ## check if geneBins already exist
  if( "geneBins" %in% names(readsOut) ){
    data.out <- c(readsOut,list('peakCallResult' = peak_call_batches) )
  }else{
    data.out <- c(readsOut,list('geneBins'=geneBins,
                                'peakCallResult' = peak_call_batches) )
  }

  return(data.out)

}


#' @title Fisher's exact test
#' @param counts.m  matrix with the first column as IP counts, and second column as input counts
#' @return data frame with p-values and odds ratio
.fisher_exact_test <- function(IP, input, IP_overall, input_overall, pseudo_count=1){

  test.m <- matrix(c(IP, input, IP_overall, input_overall), nrow = 2, byrow = FALSE,
                   dimnames = list(c("IP", "input"), c("bin", "overall")))

  # add a pseudo_count (1) to avoid zeros in denominators when calculating the odds ratio
  test.m <- test.m + pseudo_count

  fisher_test <- fisher.test(test.m, alternative="greater")

  fisher_result <- data.frame(pvalue = NA, odds_ratio = NA)

  fisher_result$pvalue <- fisher_test$p.value
  # fisher_result$odds_ratio <- fisher_test$estimate

  # use sample odds ratio (unconditional MLE) rather than the conditional Maximum Likelihood Estimate from Fisher's exact test
  a <- test.m[1,1] # IP in bin
  b <- test.m[1,2] # IP overall
  c <- test.m[2,1] # input in bin
  d <- test.m[2,2] # input overall

  # fisher_result$odds_bin <- a/c

  fisher_result$odds_ratio <- (a*d)/(b*c)

  return(fisher_result)
}
