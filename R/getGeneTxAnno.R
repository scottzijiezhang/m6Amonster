

.getGeneBins <- function(geneGRList,geneNames,binSize,threads = 1){

  registerDoParallel(cores = threads)
  geneBins <-foreach(i = 1:length(geneNames), .combine = rbind)%dopar%{
    geneName = geneNames[i]
    geneModel =reduce( geneGRList[geneName][[1]] )## merge overlapping exons

    # DNA location to gene location conversion
    df.geneModel= as.data.frame(geneModel) ##data frame of gene model
    dna.range = as.data.frame(range(geneModel) )
    df.geneModel$end = df.geneModel$end - dna.range$start + 1
    df.geneModel$start = df.geneModel$start - dna.range$start + 1
    DNA2RNA = rep(0,dna.range$end - dna.range$start +1)
    no.exon = dim(df.geneModel)[1]
    for (j in 1:no.exon){DNA2RNA[df.geneModel$start[j]:df.geneModel$end[j]]=1}
    exon.length = sum(DNA2RNA)
    #creat a corresponding map from RNA to DNA
    RNA2DNA = 1:exon.length
    pointer = 1
    for (j in 1:no.exon){
      RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1) ]= RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1)] + df.geneModel$start[j] -pointer
      pointer = pointer + df.geneModel$width[j]
    }
    RNA2DNA = RNA2DNA + dna.range$start -1 #back to chromosome coordinates
    #creat center points of continuous window
    if(exon.length <= binSize){
      slidingStart= exon.length/2
      mapping = data.frame(start = RNA2DNA[slidingStart-exon.length/2+1], end = RNA2DNA[slidingStart + exon.length/2]  )
      }else{
      slidingStart= round(seq(from = binSize/2, to = (exon.length - binSize/2), length.out = ceiling(exon.length/binSize) ) )
      mapping = data.frame(start = RNA2DNA[slidingStart - binSize/2 +1], end = RNA2DNA[slidingStart + binSize/2 ]  )
      }
    mapping$chr = as.character(dna.range$seqnames)
    mapping$strand = as.character(dna.range$strand)
    cbind(data.frame(geneName,slidingStart),mapping)
  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
  return(geneBins)
}

.getPeakBins <- function(geneGRList,geneName,slidingStarts,binSize){

    geneModel =reduce( geneGRList[geneName][[1]] )## merge overlapping exons

    # DNA location to gene location conversion
    df.geneModel= as.data.frame(geneModel) ##data frame of gene model
    dna.range = as.data.frame(range(geneModel) )
    df.geneModel$end = df.geneModel$end - dna.range$start + 1
    df.geneModel$start = df.geneModel$start - dna.range$start + 1
    DNA2RNA = rep(0,dna.range$end - dna.range$start +1)
    no.exon = dim(df.geneModel)[1]
    for (j in 1:no.exon){DNA2RNA[df.geneModel$start[j]:df.geneModel$end[j]]=1}
    exon.length = sum(DNA2RNA)
    #creat a corresponding map from RNA to DNA
    RNA2DNA = 1:exon.length
    pointer = 1
    for (j in 1:no.exon){
      RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1) ]= RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1)] + df.geneModel$start[j] -pointer
      pointer = pointer + df.geneModel$width[j]
    }
    RNA2DNA = RNA2DNA + dna.range$start -1 #back to chromosome coordinates

    #create center points of continuous window
    if(exon.length <= binSize | (dna.range$strand == "+" & slidingStarts[2] == ( exon.length - binSize - exon.length %% binSize + 1) ) ){ # for genes that is shorter than bin size || if it is the rightmost bin (buffer bin) to be retrieved (positive strand only)
      mapping = data.frame(start = RNA2DNA[ slidingStarts[1] ], end = RNA2DNA[exon.length]  )
    }else if( dna.range$strand == "-" & slidingStarts[2] == 1 ){ # when retrieve the leftmost (buffer) bin on reverse strand
      mapping = data.frame(start = RNA2DNA[ slidingStarts[1] ], end = RNA2DNA[slidingStarts[2] + binSize + exon.length %% binSize - 1 ]  )
    }else{ # retrieve regular bins
      mapping = data.frame(start = RNA2DNA[ slidingStarts[1] ], end = RNA2DNA[slidingStarts[2] + binSize - 1 ] )
    }

    mapping$chr = as.character(dna.range$seqnames)
    return(mapping[,c("chr","start","end")])

}
