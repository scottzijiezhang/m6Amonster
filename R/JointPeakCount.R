#' @title JointPeakCount
#' @param readsOut The data list as output of callPeakFisher()
#' @param joint_threshold Require a bin to be called as a peak in joint_threshold number of samples to be a joint peak.
JointPeakCount <- function(readsOut, joint_threshold = 2 ){

  if(!"peakCallResult" %in% names(readsOut)){
    stop("The readsOut input must contain peakCallResult data from callPeakFisher()...")
  }

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

  if (nrow(peak_id_pairs) > 0) {
    dist_checkpoints <- as.integer(as.character(geneBins$bin[peak_id_pairs[,"end_id"]])) - as.numeric(as.character(geneBins$bin[peak_id_pairs[,"start_id"]] ) ) # this is the original PEAK_LENGTH
    PEAK_LENGTH <- dist_checkpoints + readsOut$binSize # add a flanking window (PARAMETERS$WINDOW_WIDTH/2) to both sides of checkpoints
    # print(head(cbind(peak_id_pairs, dist_checkpoints, PEAK_LENGTH)))
    cat("Peak length distribution:\n")
    print(summary(PEAK_LENGTH))
  }

  m6A <- readsOut$reads[,(1+length(readsOut$samplenames)):(2*length(readsOut$samplenames))]
  input <- readsOut$reads[,1:length(readsOut$samplenames)]

  joint_peak_ip <- t( apply(peak_id_pairs,1,function(x,y){
    if(x[1]==x[2]){
      return(y[x[1]:x[2],])
    }else{
      colSums(y[x[1]:x[2],])
    }
  },y = m6A) )

  joint_peak_input <- t( apply(peak_id_pairs,1,function(x,y){
    if(x[1]==x[2]){
      return(y[x[1]:x[2],])
    }else{
      colSums(y[x[1]:x[2],])
    }
  },y = input) )

  data.out <- c(readsOut,
                list('jointPeak_id_pairs' = peak_id_pairs,
                     'jointPeak_ip' = joint_peak_ip,
                     'jointPeak_input' = joint_peak_input)
                )
  return(data.out)

}
