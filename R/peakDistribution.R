#' @title gtfToGeneModel
#' @description to remove ambiguous gene model and return gene model as genomic ranges object
#' @param bedFile peak list as bed format
#' @param gtf gtf file to build gene model
#' @param saveName the file name to save ditribution plot
#' @export
peakDistribution <- function(peak,gtf,saveName = NA){

  library(GenomicFeatures)


  ## collapes the peak to the center
  x <- peak
  for(i in 1:dim(x)[1]){
    start = x[i,2]
    end = x[i,3]
    blocks = x[i,10]
    blockStart = as.numeric( unlist(strsplit(as.character(x[i,12]),",")) )
    blockSize = as.numeric( unlist(strsplit(as.character(x[i,11]),",")) )
    if(blocks <2){
      x[i,2] =  x[i,3] = round(start + sum(blockSize)/2 )
    } else{
      pointer = 2
      while(sum(blockSize[1:pointer]) < sum(blockSize)/2 ){pointer = pointer+1}
      x[i,2] =  x[i,3] = round(start + blockStart[pointer] + sum(blockSize)/2 - sum(blockSize[1:(pointer-1)]) )
    }
    x[i,10] = 1
  }

  ################
  gr.peak <- makeGRangesFromDataFrame(df = x,start.field = "V2",end.field = "V3",seqnames.field = "V1",strand.field = "V6")
  txdb <- makeTxDbFromGFF(file = gtf,format = "gtf")

  cds <-  cdsBy(txdb,by = "tx")
  n.cds <- length( which(countOverlaps(gr.peak,cds)>0) )

  fiveUTR <- fiveUTRsByTranscript(txdb)
  n.fiveUTR <- length( which(countOverlaps(gr.peak,fiveUTR)>0) )

  threeUTR <- threeUTRsByTranscript(txdb)
  n.threeUTR <-  length( which(countOverlaps(gr.peak,threeUTR)>0) )

  exon <- exonsBy(txdb,by = "tx")
  all_mRNA <- unique(c(names(fiveUTR),names(threeUTR),names(cds)))
  name_ncRNA <- setdiff(names(exon),all_mRNA)
  ncRNA <- exon[name_ncRNA]
  n.ncRNA <-  length( which(countOverlaps(gr.peak,ncRNA)>0) )

  slices <- c(n.cds,n.fiveUTR,n.threeUTR,n.ncRNA)
  lbls <- c("CDS","5'UTR","3'UTR", "ncRNA")
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(lbls, pct) # add percents to labels
  lbls <- paste(lbls,"%",sep="") # ad % to labels
  if(is.na(saveName)){
    pie(slices,labels = lbls, col= rainbow(4), main="Pie Chart of m6A peak distribution")
  } else{
    pdf(saveName)
    pie(slices,labels = lbls, col=rainbow(4), main="Pie Chart of m6A peak distribution")
    dev.off()
  }

}
