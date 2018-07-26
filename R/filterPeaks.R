#' @title filterPeak
#' @param readsOut The data list contain peak calling result
#' @param minINPUT Define the minimal number of INPUT reads (median across samples).
#' @param minIP Define the minimal number of IP reads (median across samples)
#' @return readsOut with updated jointPeak_id_pairs that have been filtered by read count.
#' @export
filterPeak <- function(readsOut, minINPUT = 15, minIP = 30){
  ## check if joint peak called
  if("norm.jointPeak_ip"%in%names(readsOut) & "jointPeak_id_pairs"%in%names(readsOut) & "jointPeak_ip"%in%names(readsOut)& "jointPeak_input"%in%names(readsOut)){
    input <- readsOut$jointPeak_input
    ip <- readsOut$jointPeak_ip

    keep.id <- intersect(  which(rowMedians(input) >= minINPUT ), which(rowMedians(ip) >= minIP ) )
    keep.id <- sort(keep.id)

    data.out <- c(readsOut[1:8],
                  list('jointPeak_id_pairs' = readsOut$jointPeak_id_pairs[keep.id, ],
                       'jointPeak_ip' = ip[keep.id, ],
                       'jointPeak_input' = input[keep.id, ]),
                  readsOut["geneSum"],
                  list('norm.jointPeak_ip'= readsOut$norm.jointPeak_ip[keep.id, ]),
                  readsOut[c("sizeFactor","X")]
    )

    if("jointPeak_adjExpr"%in%names(readsOut) ){
      cat("Please run adjustExprLevel() again after filtering peaks...")
    }

    return(data.out)

  }else if("jointPeak_id_pairs"%in%names(readsOut) & "jointPeak_ip"%in%names(readsOut)& "jointPeak_input"%in%names(readsOut)){
    input <- readsOut$jointPeak_input
    ip <- readsOut$jointPeak_ip

    keep.id <- unique( c( which(rowMedians(input) >= minINPUT ), which(rowMedians(ip) >= minIP ) ) )
    keep.id <- sort(keep.id)

    data.out <- c(readsOut[1:8],
                  list('jointPeak_id_pairs' = readsOut$jointPeak_id_pairs[keep.id, ],
                       'jointPeak_ip' = ip[keep.id, ],
                       'jointPeak_input' = input[keep.id, ])
                  )
    return(data.out)

  }else{
    stop("Needs to get joint peak read count before calling filterPeak()...")
  }



}
