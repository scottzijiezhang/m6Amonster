######################
# The main function ##
######################
peakQuant <-
  function(
    bedFiles,
    filenames,# file name of samples
    gtf, # gtf file used for peak calling
    shift = 75, #number of nucleotide shifted at peak calling
    modification = "m6A",
    sharedInput=FALSE,
    inputNames=NA
  ){
    if( !all(file.exists(bedFiles)) ) stop( "one or all bed files missing" )
    ## make the tx database object
    txdb=makeTxDbFromGFF(gtf,format="gtf")
    print("obtained gene model coordinates...")

    dir.create(paste(modification,"_peak_quantification",sep = ""))
    print(paste( "Creating folder ",getwd(),"/" ,modification,"_peak_quantification/",sep = ""))

    mergedPeak = mergePeak(bedFiles = bedFiles, gtf = txdb , outputDir = paste(modification,"_peak_quantification/",sep = ""))

    countPeakReads(filenames = filenames,mergedPeak = mergedPeak, gtf = txdb, shift = shift, modification=modification, sharedInput=sharedInput, inputNames=inputNames )

    print("peak quantification... done done done.. done.....")

}

#######################################
## a helper function to merge peaks ###
#######################################
mergePeak <-
  function(bedFiles,gtf,outputDir){
    ## make the tx database object
    if(is.character(gtf)){
      txdb=makeTxDbFromGFF(gtf,format="gtf")
    }else{
      txdb = gtf
    }
    exons=exonsBy(txdb,by="gene")
    print("obtained gene model coordinates...")

    ## bind all peaks into one data frame.
    peaks = data.frame()
    for (i in 1:length(bedFiles)){
      tmp=read.table(bedFiles[i],header=F)
      peaks=rbind(peaks,tmp)
    }
    peaks=peaks[,c(1,2,3,6,4)]## retain bed4 info
    colnames(peaks)=c("chr","start","end","strand","names")
    ##convert the dataframe into genomic ranges
    peak.gr=makeGRangesFromDataFrame(peaks,keep.extra.columns=T)
    ## remove redundancy by merging overlapping peaks
    peak.gr.merge=reduce(peak.gr)
    map=as.matrix( findOverlaps(peak.gr.merge,peak.gr) )
    ID=map[!duplicated(map[,1]),][,2]#ID of reduced peak in original peak.gr
    ## retrieve the gene symbol of the peak
    peak.gr.merge$names=peak.gr[ID]$names
    print("all overlapping peak merged...")

    peak.genes=exons[unique(peak.gr.merge$names)]##subsetting the annotation for peaked genes
    no.peaks=length(peak.gr.merge)## define number of total peaks
    bed12=data.frame(matrix(nrow=no.peaks,ncol=12)) # initiate the bed12 data frame to store the peaks
    colnames(bed12)=c("chr","start","end","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSize","blockStart")

    print("preparing bed12 files for the merged peak...")
    for(i in 1:no.peaks){
      tmp=peak.gr.merge[i]
      geneModel=peak.genes[tmp$names]
      tmp=GenomicRanges::intersect(tmp,peak.genes[tmp$names][[1]])
      tmp$names = names(geneModel)
      bed12[i,"chr"] = unique(as.character(seqnames(tmp)))
      bed12[i,c(2,3)] = as.data.frame(range(tmp))[1,c(2,3)]
      bed12[i,"strand"]= unique(as.character(strand(tmp)))
      bed12[i,"name"] = names(geneModel)
      bed12[i,"score"] = 1
      bed12[i,c("thickStart","thickEnd")] = bed12[i,c(2,3)]
      bed12[i,"itemRgb"] = NA
      bed12[i,"blockCount"] = length(tmp)
      bed12[i,"blockSize"] = paste(as.data.frame(tmp)[,4],collapse=",")
      bed12[i,"blockStart"] = paste(as.data.frame(tmp)[,2] - replicate(bed12[i,"blockCount"],bed12[i,2]),collapse=",")
    }

    ##write the bed12 file to the specified directory
    write.table(bed12,file=paste(outputDir,"mergedPeak.bed",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

    return(peak.gr.merge)

}

######################################
## A helper function to count reads ##
######################################
countPeakReads <-
  function(
    filenames,# file name of samples
    mergedPeak,# merged peak as genomic range object
    gtf, # gtf file used for peak calling
    shift = 75, #number of nucleotide shifted at peak calling
    modification = "m6A",
    sharedInput=FALSE,
    inputNames=NA
  ){
    ##read bam files
    if(sharedInput){
      bamPath.input = paste(inputNames,".input.bam",sep="")
      bamPath.IP = paste(filenames,".",modification,".bam",sep="")
      no.samples = length(filenames)
    }else{
      bamPath.input = paste(filenames,".input.bam",sep="")
      bamPath.IP = paste(filenames,".",modification,".bam",sep="")
      no.samples = length(filenames)
    }


    if( !all(file.exists(bamPath.input)) ) stop( "input bam file missing!!!" )
    if( !all(file.exists(bamPath.IP)) ) stop( "IP bam file missing!!!" )

    #load gene models for peak genes
    ## make the tx database object
    if(is.character(gtf)){
      txdb=makeTxDbFromGFF(gtf,format="gtf")
    }else{
      txdb = gtf
    }
    exons=exonsBy(txdb,by="gene")
    peak.genes=exons[unique(mergedPeak$names)]


    no.peaks=length(mergedPeak)## define number of total peaks
    peakCount.input=as.data.frame( matrix(nrow=no.peaks,ncol=no.samples + 3) ) #initiate the matrix to store read count in peak window of input
    peakCount.IP=as.data.frame( matrix(nrow=no.peaks,ncol=no.samples + 3 ) )#initiate the matrix to store read count in peak window of IP
    medCount.input=as.data.frame( matrix(nrow=no.peaks,ncol=no.samples + 3 ) )#initiate the matrix to store read count in none-peak window of input
    medCount.IP=as.data.frame( matrix(nrow=no.peaks,ncol=no.samples + 3 ) ) #initiate the matrix to store read count in none-peak window of IP
    colnames(peakCount.input)=colnames(peakCount.IP)=colnames(medCount.input)=colnames(medCount.IP)=c("gene","startOnGene","endOnGene",filenames)

    print("counting reads for each peaks, this step may takes a few hours....")
    for(i in 1:no.peaks){
      tmp=mergedPeak[i]
      gene.name = as.character(tmp$names)
      ## get the longest isoform gene model for the peak
      geneModel= reduce(peak.genes[ gene.name ][[1]])
      geneModel= geneModel[strand(geneModel)==as.character(strand(tmp[1])) ]## exclude exons on different strand
      geneModel= geneModel[seqnames(geneModel)==as.character(seqnames(tmp[1])) ] ## exlude exons on different chromosome
      ## drop unused seqlevels
      seqlevels(geneModel,force=T)=seqlevels(mergedPeak)

      # DNA location to gene location conversion
      df.geneModel= as.data.frame(geneModel) ##data frame of gene model
      dna.range = as.data.frame(range(geneModel))
      df.geneModel$end = df.geneModel$end - dna.range$start + 1
      df.geneModel$start = df.geneModel$start - dna.range$start + 1
      DNA2RNA = rep(0,dna.range$end - dna.range$start +1)
      no.exon = dim(df.geneModel)[1]
      for (j in 1:no.exon){DNA2RNA[df.geneModel$start[j]:df.geneModel$end[j]]=1}
      exon.length = sum(DNA2RNA)
      DNA2RNA=cumsum(DNA2RNA)*DNA2RNA

      ## obtain peak in RNA coordinate
      peak = as.data.frame(tmp)[,c("start","end","strand")]
      if(peak$start <= dna.range$start ){ peak$start = dna.range$start } #check out of boundary peak
      peak$start = DNA2RNA[peak$start-dna.range$start + 1 ]
      if (peak$end > dna.range$end  ){ peak$end = dna.range$end } #check out of boundary peak
      peak$end = DNA2RNA[peak$end-dna.range$start + 1 ]
      if(peak$strand == "+"){peak$strand = "-"}else{peak$strand = "+"} ## switch strand on RNA

      #sliding window
      if(exon.length<=100){
        slidingStart= 0
      }else{
        slidingStart= round(seq(from = 50, to = (exon.length-50), length.out = (exon.length-100)/10 ) )
      }


      #count reads in all samples
      ba.IP = sapply(bamPath.IP,countReadFromBam,which = range(geneModel),peak=peak,DNA2RNA = DNA2RNA,shift=shift,left=dna.range$start,sliding = slidingStart)
      ba.input = sapply(bamPath.input,countReadFromBam,which = range(geneModel),peak=peak,DNA2RNA = DNA2RNA,shift=shift,left=dna.range$start,sliding = slidingStart)

      #write the result to global variables
      peakCount.input[i,"gene"]=peakCount.IP[i,"gene"]=medCount.input[i,"gene"]=medCount.IP[i,"gene"] = gene.name
      peakCount.input[i,"startOnGene"]=peakCount.IP[i,"startOnGene"]=medCount.input[i,"startOnGene"]=medCount.IP[i,"startOnGene"] = peak$start
      peakCount.input[i,"endOnGene"]=peakCount.IP[i,"endOnGene"]=medCount.input[i,"endOnGene"]=medCount.IP[i,"endOnGene"] = peak$end
      peakCount.IP[i,4:(no.samples+3)] = ba.IP[1,]
      peakCount.input[i,4:(no.samples+3)] = ba.input[1,]
      medCount.IP[i,4:(no.samples+3)] = ba.IP[2,]
      medCount.input[i,4:(no.samples+3)] = ba.input[2,]
      print(i)
    }

    write.table(peakCount.input,file = paste(modification, "_peak_quantification/", modification,"_peakCount.input.xls",sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)
    write.table(medCount.input,file = paste(modification, "_peak_quantification/", modification,"_medCount.input.xls",sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)
    write.table(peakCount.IP,file = paste(modification, "_peak_quantification/", modification,"_peakCount.IP.xls", sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)
    write.table(medCount.IP,file = paste(modification, "_peak_quantification/", modification,"_medCount.IP.xls", sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)
    save(peak.genes,peakCount.input,peakCount.IP,medCount.input,medCount.IP,mergedPeak,file = paste(modification, "_peak_quantification/peakQuant.RData",sep = "") )
    print(paste("write counts files to ", getwd(),"/",modification, "_peak_quantification/","...", sep = "" ))
  }

##################################
## A helper function of helper ###
##################################
countReadFromBam <-
  function(bam,which,peak,DNA2RNA,shift,left,sliding){
    ba = scanBam(bam, param=ScanBamParam( which=which, what =c("pos","strand") ) )
    ba = data.frame( pos = ba[[1]]$pos, strand = ba[[1]]$strand )
    ba = ba[ ba$pos > left, ]
    ba = ba[ba$strand == peak$strand, ] ## filter for strand
    ba$pos = DNA2RNA[ba$pos - left] ## convert mapped read pos into RNA position
    ba = ba[ ba$pos >0, ] ## drop intron reads.
    ##shift the read pos to the center of the reads
    if(peak$strand == "+"){ba$pos = ba$pos + shift}else{ba$pos = ba$pos +50 - shift}
    ##count the reads in the sliding windows
    no.window = length(sliding)
    windowCounts = vector(length = no.window)
    for(j in 1:no.window){
      windowCounts[j]= sum( abs(sliding[j] - ba$pos) < 50 )
    }
    medCounts = median(windowCounts) # median coounts of all sliding window on this gene

    peakCounts = sum( ba$pos >= peak$start & ba$pos <= peak$end )
    return(c(peakCounts,medCounts))
}



