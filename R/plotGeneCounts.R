#' @title plotGeneCounts
#' @param IP_BAM
#' @param INPUT_BAM
#' @param size.IP
#' @param size.INPUT
#' @param geneName
#' @param geneModel
#' @param libraryType can be "opposite" or "same"
## the main function to plot m6A-seq on one group of data
plotGeneCounts <- function(IP_BAM, INPUT_BAM, size.IP, size.INPUT, geneName, geneModel, libraryType = "opposite", center = mean ,ZoomIn=NULL){
  IP.cov <- getAveCoverage(geneModel= geneModel,bamFiles = IP_BAM,geneName = geneName,size.factor = size.IP,center = center, ZoomIn = ZoomIn)
  INPUT.cov <- getAveCoverage(geneModel= geneModel,bamFiles = INPUT_BAM,geneName = geneName,size.factor = size.INPUT, center = center,ZoomIn = ZoomIn)
  cov.data <- data.frame(IP=IP.cov,Input=INPUT.cov,genome_location=as.numeric(names(IP.cov) ) )
  ggplot(data = cov.data,aes(genome_location))+geom_line(aes(y=Input,colour ="Input"))+geom_line(aes(y=IP,colour="IP"))+labs(y="normalized coverage")+scale_x_continuous(breaks = round(seq(min(cov.data$genome_location), max(cov.data$genome_location), by = ((max(cov.data$genome_location)-min(cov.data$genome_location))/10) ),1))
}


#' @title plotGenePairCounts
#' @description plot tow groups of samples in the same figure
#' @param Ctl_IP_BAM
#' @param Ctl_INPUT_BAM
#' @param Treat_IP_BAM
#' @param Treat_INPUT_BAM
#' @param Ctl_size.IP
#' @param Ctl_size.INPUT
#' @param Treat_size.IP
#' @param Treat_size.INPUT
#' @param geneName
#' @param geneModel
## the main function to plot m6A-seq on two group of data
plotGenePairCounts <- function(Ctl_IP_BAM,Ctl_INPUT_BAM,Treat_IP_BAM,Treat_INPUT_BAM,Ctl_size.IP,Ctl_size.INPUT,Treat_size.IP,Treat_size.INPUT,geneName,geneModel, libraryType = "opposite",center = mean,ZoomIn=NULL){
  Ctl_IP.cov <- getAveCoverage(geneModel= geneModel,bamFiles = Ctl_IP_BAM,geneName = geneName,size.factor = Ctl_size.IP,center = center, ZoomIn = ZoomIn)
  Ctl_INPUT.cov <- getAveCoverage(geneModel= geneModel,bamFiles = Ctl_INPUT_BAM,geneName = geneName,size.factor = Ctl_size.INPUT,center = center , ZoomIn = ZoomIn)
  Treat_IP.cov <- getAveCoverage(geneModel= geneModel,bamFiles = Treat_IP_BAM,geneName = geneName,size.factor = Treat_size.IP, center = center,ZoomIn = ZoomIn)
  Treat_INPUT.cov <- getAveCoverage(geneModel= geneModel,bamFiles = Treat_INPUT_BAM,geneName = geneName,size.factor = Treat_size.INPUT, center = center,ZoomIn = ZoomIn)
  cov.data <- data.frame(Ctl_IP=Ctl_IP.cov, Ctl_Input = Ctl_INPUT.cov,
                         Treat_IP=Treat_IP.cov, Treat_Input = Treat_INPUT.cov,
                         genome_location=as.numeric(names(Ctl_IP.cov) ) )
  ggplot(data = cov.data,aes(genome_location))+geom_line(aes(y=Ctl_Input,colour ="Ctl Input"))+geom_line(aes(y=Treat_IP,colour="Treat IP"))+geom_line(aes(y=Treat_Input,colour ="Treat Input"))+geom_line(aes(y=Ctl_IP,colour="Ctl IP"))+labs(y="normalized coverage")+scale_x_continuous(breaks = round(seq(min(cov.data$genome_location), max(cov.data$genome_location), by = ((max(cov.data$genome_location)-min(cov.data$genome_location))/10) ),1))
}

## helper function to get average coverage of a gene of multiple samples
getAveCounts <- function(geneModel,bamFiles,geneName,size.factor, center , libraryType = "opposite" ,ZoomIn){
  locus <- as.data.frame( range(geneModel[geneName][[1]]) )
  if(is.null(ZoomIn)){
  }else{
    locus$start = ZoomIn[1]
    locus$end = ZoomIn[2]
    locus$width = ZoomIn[2] - ZoomIn[1] + 1
  }
  covs <- sapply(bamFiles,getCounts,locus=locus, libraryType = libraryType)
  covs <- covs/size.factor
  ave.cov <- apply(covs,1, center)
  return(ave.cov)
}

getCounts <- function(bf,locus,libraryType){
  s_param <- ScanBamParam(which = GRanges(locus$seqnames,IRanges(locus$start,locus$end)),what = c("strand", "pos"), mapqFilter = 30 )
  #p_param <- PileupParam(max_depth=1000000,min_nucleotide_depth=0,distinguish_nucleotides=F,)
  #get coverage from the bam file
  res <- scanBam(bf,param = s_param)
  res <- as.data.frame(res)
  if(libraryType == "opposite"){
    res <- res[res[,1]!=locus$strand,]
  }else if (libraryType == "same"){
    res <- res[res[,1]==locus$strand,]
  }else{
    stop("libraryType must be opposite or same... ")
  }
  counts <- table(res[,2])
  cov <- vector(length = locus$width)
  names(cov) <- c(locus$start:locus$end)
  cov[1:locus$width] <- 0
  cov[names(counts)] <- counts
  return(cov)
}
