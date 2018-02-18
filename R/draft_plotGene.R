plotGene <- function(IP_BAM, INPUT_BAM, size.IP, size.INPUT, geneName, geneModel, libraryType = "opposite", center = mean ,ZoomIn=NULL,gtf){
  IP.cov <- getAveCoverage(geneModel= geneModel,bamFiles = IP_BAM,geneName = geneName,size.factor = size.IP, libraryType = libraryType, center = center, ZoomIn = ZoomIn)
  INPUT.cov <- getAveCoverage(geneModel= geneModel,bamFiles = INPUT_BAM,geneName = geneName,size.factor = size.INPUT, libraryType = libraryType, center = center,ZoomIn = ZoomIn)
  cov.data <- data.frame(IP=IP.cov,Input=INPUT.cov,genome_location=as.numeric(names(IP.cov) ) )
  yscale <- max(IP.cov,INPUT.cov)
  p1 <- "ggplot(data = cov.data,aes(genome_location))+geom_line(aes(y=Input,colour =\"Input\"))+geom_line(aes(y=IP,colour=\"IP\"))+labs(y=\"normalized coverage\")+scale_x_continuous(breaks = round(seq(min(cov.data$genome_location), max(cov.data$genome_location), by = ((max(cov.data$genome_location)-min(cov.data$genome_location))/10) ),1))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = \"black\"))"

  GTF <- rtracklayer::import(gtf,format = "gtf")
  cds <- GTF[GTF$type == "CDS" & GTF$gene_name == geneName]
  p2 <- .getGeneModelAnno(geneModel,geneName,cds)
  p <- paste(p1,p2,sep = "+")
  eval(parse( text = p ))
}
