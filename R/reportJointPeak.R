#' @title reportJointPeak
#' @param readsOut The data list contain peak calling result
#' @param joint_threshold Define the number of sample required to have consistent peak in a locus to call joint peak.
#'  This option is useful only when joint poeak has not been reported by jointPeakCount().
#' @return merged.report The joint peak of merged bins.
reportJointPeak <- function(readsOut, joint_threshold,threads = 1){
  if("all.est.peak" %in% names(readsOut) & "peakCallResult"%in% names(readsOut)){
    cat("Reporting joint peak with PoissonGamma test statistics...\n")
    peak_id_pairs <- readsOut$jointPeak_id_pairs

    stats <- readsOut$all.est.peak
    geneGRList <- readsOut$geneModel
    peakGenes <- as.character(readsOut$geneBins[peak_id_pairs[,1],"gene"])
    ##
    if("p_value3" %in% colnames(stats)){ colnames(stats)[which(colnames(stats) == "p_value3")] = "p_value"}

    num_peaks <- nrow(peak_id_pairs)
    cat(paste("Reporting ",num_peaks," peaks...\n"))
    if (num_peaks == 0){return(data.frame())
    }else {
      start_time <- Sys.time()
      registerDoParallel(cores = threads)
      cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
      cat(paste("Using",getDoParWorkers(),"thread(s) to report merged report...\n"))
      merged.report<- foreach( p = 1:num_peaks, .combine = rbind)%dopar%{
        peak_row_id <- peak_id_pairs[p,]
        geneExons <- geneGRList[peakGenes[p]][[1]]

        peak <- .getPeakBins(geneGRList,peakGenes[p],c(readsOut$geneBins$bin[peak_row_id[1]],readsOut$geneBins$bin[peak_row_id[2]]),readsOut$binSize )
        peakE <- .peakExons(peak,as.data.frame(geneExons))
        data.frame(chr=peak$chr,
                   start = peak$start,
                   end = peak$end,
                   name = peakGenes[p],
                   score = 0,
                   strand = as.character(strand(geneExons))[1],
                   thickStart = peak$start,
                   thickEnd = peak$end,
                   itemRgb=0,
                   blockCount = nrow(peakE),
                   blockSizes = paste(peakE$width,collapse=","),
                   blockStarts = paste(peakE$start - replicate(nrow(peakE),peakE$start[1]),collapse=","),
                   logFC = stats[p,"beta"],
                   p_value = stats[p, "p_value"],
                   fdr = stats[p,"padj"]
        )
      }
      rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
      end_time <- Sys.time()
      cat(paste("Time used to report peaks:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
    }
    return( merged.report )

  }else if("peakCallResult"%in% names(readsOut)){
    cat("Reporting joint peak without differential peak test...\n")
    ## Get joint peak
    geneBins <- readsOut$geneBins
    ## set logic vector for joint peak
    ID <- (rowSums(readsOut$peakCallResult) >= joint_threshold)
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

    peak_id_pairs <- cbind(start_id, end_id)

    geneGRList <- readsOut$geneModel
    peakGenes <- as.character(geneBins[peak_id_pairs[,1],"gene"])
    #geneBins <- .getGeneBins(geneGRList,peakGenes,readsOut$binSize )

    if (num_peaks == 0){return(data.frame())
    }else {
      start_time <- Sys.time()
      registerDoParallel(cores = 4)
      cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
      cat(paste("Using",getDoParWorkers(),"thread(s) to report merged report...\n"))
      merged.report<- foreach( p = 1:num_peaks, .combine = rbind)%dopar%{
        peak_row_id <- peak_id_pairs[p,]
        geneExons <- geneGRList[peakGenes[p]][[1]]

        peak <- .getPeakBins(geneGRList,peakGenes[p],c(geneBins$bin[peak_row_id[1]],geneBins$bin[peak_row_id[2]]),readsOut$binSize )

        peakE <- .peakExons(peak,as.data.frame(geneExons))
        data.frame(chr=peak$chr,
                   start = peak$start,
                   end = peak$end,
                   name = peakGenes[p],
                   score = 0,
                   strand = as.character(strand(geneExons))[1],
                   thickStart = peak$start,
                   thickEnd = peak$end,
                   itemRgb=0,
                   blockCount = nrow(peakE),
                   blockSizes = paste(peakE$width,collapse=","),
                   blockStarts = paste(peakE$start - replicate(nrow(peakE),peakE$start[1]),collapse=",")
        )
      }
      rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
      end_time <- Sys.time()
      cat(paste("Time used to report peaks:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
    }
    return( merged.report )

  }


}


.peakExons <- function(peak,y){
  exonID <- peak$start <= y$end & peak$end >= y$start
  if(sum(exonID) == 1){
    return(data.frame(start = peak$start, end = peak$end, width = peak$end - peak$start))
  }else if(sum(exonID) > 1){
    peakexon <- y[exonID,]
    peakexon[1,"start"] <- peak$start
    peakexon[sum(exonID),"end"] <- peak$end
    return(data.frame(start = peakexon$start, end = peakexon$end, width = peakexon$end - peakexon$start))
  }
}
