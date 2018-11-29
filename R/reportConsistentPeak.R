#' @title reportConsistentPeak
#' @param readsOut The data list contain peak calling result
#' @param samplenames The samplenames to be reported for consistent peaks
#' @param joint_threshold Define the number of sample required to have consistent peak in a locus to call consistent peak in a group.
#'  This option is useful only when joint poeak has not been reported by jointPeakCount().
#' @return merged.report The joint peak of merged bins.
#' @export
reportConsistentPeak <- function(readsOut, samplenames ,joint_threshold=NULL,threads = 1){
  if(! all(samplenames %in% readsOut$samplenames ) ){
    stop("Not all samples in dataset!")
  }else{
    ## define samples to be reported
    sample_report_id <- match(samplenames, readsOut$samplenames )

    if(is.null(joint_threshold)){
      joint_threshold <- length(samplenames)
      cat("Reporting peak concsistent in all samples for\n",paste(samplenames,collapse = " "),"\n")
    }else if(joint_threshold > length(samplenames) ){
      joint_threshold <- length(samplenames)
      cat("Reporting peak concsistent in all samples for\n",paste(samplenames,collapse = " "),"\n")
    }else{
      cat("Reporting joint peak consistant in at least ",joint_threshold," samples among \n",paste(samplenames,collapse = " "),"\n")
    }
  }

  if("peakCallResult" %in% names(readsOut)){

    ## Get joint peak
    geneBins <- readsOut$geneBins

    ## set logic vector for joint peak
    if(length(sample_report_id) >1 ){
      ID <- (rowSums(readsOut$peakCallResult[,sample_report_id]) >= joint_threshold)
    }else if(length(sample_report_id) == 1){
      ID <- readsOut$peakCallResult[,sample_report_id]
    }

    num_lines <- length(ID)

    # start ids of checkpoints
    ## find peak-starting checkpoints (either from nonpeak to peak, or peak in a different batch)
    start_id <- which((ID[2:num_lines]-ID[1:num_lines-1]==1) |
                        ((geneBins$gene[2:num_lines]!=geneBins$gene[1:num_lines-1]) & (ID[2:num_lines] == TRUE)) )
    start_id <- start_id + 1 # add 1 since ID was counted from 2 to num_lines
    if ( ID[1]==TRUE ) { start_id <- c(1,start_id) } # if the first checkpoint bin is peak

    # end ids of checkpoints
    ## find peak-ending checkpoints (either from peak to nonpeak, or peak in a different batch)
    end_id <- which((ID[1:num_lines-1]-ID[2:num_lines]==1) |
                      ((geneBins$gene[1:num_lines-1]!=geneBins$gene[2:num_lines]) & (ID[1:num_lines-1] == TRUE)) )
    if (ID[num_lines]==TRUE) {end_id <- c(end_id,num_lines)} # if the last checkpoint bin is peak

    ## get peak id-pairs that defines the bins to be merged as one peak
    peak_id_pairs <- cbind(start_id, end_id)
    num_peaks <- nrow(peak_id_pairs)
    geneGRList <- readsOut$geneModel
    peakGenes <- as.character(geneBins[peak_id_pairs[,1],"gene"])

    ## Get raw readcount for peaks.
    m6A <- readsOut$reads[,(length(readsOut$samplenames)+sample_report_id)]
    input <- readsOut$reads[,sample_report_id]
    ## standardize library size
    all.size <- colSums( cbind(input,m6A) )/mean( colSums( cbind(input,m6A) ) )

    if(length(sample_report_id) >1 ){
      ## get Summed read count
      ip_sum <- round( rowSums( t( t(m6A)/all.size[ length(sample_report_id)+1:length(sample_report_id) ] ) ) )
      input_sum <- round(rowSums( t( t(input)/all.size[ 1:length(sample_report_id) ] ) ) )
    }else if(length(sample_report_id) == 1){
      ## No need to sum for single sample
      ip_sum <- round( m6A/all.size[2] )
      input_sum <- round( input/all.size[1] )
    }


    ## Get peak ip count
    joint_peak_ip <- apply(peak_id_pairs,1,function(x,y){
      if(x[1]==x[2]){
        return(y[x[1]:x[2]])
      }else{
        sum(y[x[1]:x[2]])
      }
    },y = ip_sum)

    ## Get peak input count
    joint_peak_input <- apply(peak_id_pairs,1,function(x,y){
      if(x[1]==x[2]){
        return(y[x[1]:x[2]])
      }else{
        sum(y[x[1]:x[2]])
      }
    },y = input_sum)

    ## Get peak gene ip median
#    joint_peak_ip_median <-  round( tapply(ip_sum,geneBins$gene,median)[peakGenes]  )
    ## Get peak gene input median
#    joint_peak_input_median <- round( tapply(input_sum,geneBins$gene,median)[peakGenes]  )

    ## Compute fisher's test for each peak
#    peak_test <- do.call( rbind.data.frame, apply(cbind(joint_peak_ip,joint_peak_input,joint_peak_ip_median,joint_peak_input_median), 1, function(x){.fisher_exact_test(x[1],x[2],x[3],x[4])}) )
    peak_test <- data.frame(pvalue = .Bino_test( joint_peak_ip, joint_peak_input, sum(ip_sum), sum(input_sum) ),
                            odds_ratio = (joint_peak_ip/sum(ip_sum))/(joint_peak_input/sum(input_sum)) )

    if (num_peaks == 0){return(data.frame())
    }else {
      start_time <- Sys.time()
      registerDoParallel(cores = threads)
      cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
      cat(paste("Using",getDoParWorkers(),"thread(s) to report merged report...\n"))
      merged.report<- foreach( p = 1:num_peaks, .combine = rbind)%dopar%{
        peak_row_id <- peak_id_pairs[p,]
        geneExons <- reduce ( geneGRList[peakGenes[p]][[1]] )

        peak <- .getPeakBins(geneGRList,peakGenes[p],c(geneBins$bin[peak_row_id[1]],geneBins$bin[peak_row_id[2]]),readsOut$binSize )

        peakE <- .peakExons(peak,as.data.frame(geneExons))
        data.frame(chr=peak$chr,
                   start = peak$start,
                   end = peak$end,
                   name = peakGenes[p],
                   score = peak_test$pvalue[p],
                   strand = as.character(strand(geneExons))[1],
                   thickStart = peak$start,
                   thickEnd = peak$end,
                   itemRgb=0,
                   blockCount = nrow(peakE),
                   blockSizes = paste(peakE$width,collapse=","),
                   blockStarts = paste(peakE$start - replicate(nrow(peakE),peakE$start[1]),collapse=","),
                   enrichmentScore = peak_test$odds_ratio[p]
        )
      }
      rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
      end_time <- Sys.time()
      cat(paste("Time used to report peaks:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
    }

    return( merged.report )

  } else{
    stop("Please run callPeakFisher() before reporting peaks...")
  }

}
