
#' @title plotExpression
#' @param RDM The RDM object
#' @param geneName The name of genes to be ploted.
#' @param logCount where to plot count at log scale
#' @export
plotExpression <- function(RDM, geneName, logCount = FALSE){
  if(length(geneName) ==1){
    temp <- as.data.frame(t(RDM$geneSum[geneName,] ) )
  }else{
    temp <- as.data.frame(RDM$geneSum[geneName,] )
  }

  if(logCount){
    temp <- log(temp)
    colnames(temp) <- paste0(RDM$X,1:length(RDM$X))
    temp$name <- factor(geneName,levels = geneName)
    temp_melt <- reshape2::melt(temp,id.vars = "name")
    temp_melt$Group <- unique(RDM$X)[1]
    for(i in 2:length(RDM$X)){
      temp_melt$Group[grep(unique(RDM$X)[i],temp_melt$variable)] <- unique(RDM$X)[i]
    }

    axis.font <- element_text(face = "bold", color = "black")
    ggplot(temp_melt, aes(x=name,y=value,fill=Group))+geom_boxplot()+labs(x="Gene Symbol",y="Log normalized read counts")+
      theme(axis.title =axis.font, axis.text = axis.font)+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         axis.line = element_line(colour = "black",size = 1),
                         axis.title.x=element_text(size=20, face="bold", hjust=0.5,family = "arial"),
                         axis.title.y=element_text(size=20, face="bold", vjust=0.5, angle=90,family = "arial"),
                         legend.title=element_text(size = 15,face = "bold"),legend.text = element_text(size = 18, face = "bold",family = "arial"),
                         axis.text.x = element_text(size = 15,face = "bold",family = "arial",colour = "black") ,axis.text.y = element_text(size = 15,face = "bold",family = "arial"),
                         plot.title = element_text(size=22, face="bold", hjust=0.5,vjust=0.5,family = "arial"))+
      ggtitle("Gene expression level")
  }else{
    colnames(temp) <- paste0(RDM$X,1:length(RDM$X))
    temp$name <- factor(geneName,levels = geneName)
    temp_melt <- reshape2::melt(temp,id.vars = "name")
    temp_melt$Group <- unique(RDM$X)[1]
    for(i in 2:length(RDM$X)){
      temp_melt$Group[grep(unique(RDM$X)[i],temp_melt$variable)] <- unique(RDM$X)[i]
    }
    axis.font <- element_text(face = "bold", color = "black")
    ggplot(temp_melt, aes(x=name,y=value,fill=Group))+geom_boxplot()+labs(x="Gene Symbol",y="Normalized read counts")+
      theme(axis.title =axis.font, axis.text = axis.font)+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         axis.line = element_line(colour = "black",size = 1),
                         axis.title.x=element_text(size=20, face="bold", hjust=0.5,family = "arial"),
                         axis.title.y=element_text(size=20, face="bold", vjust=0.5, angle=90,family = "arial"),
                         legend.title=element_text(size = 15,face = "bold"),legend.text = element_text(size = 18, face = "bold",family = "arial"),
                         axis.text.x = element_text(size = 15,face = "bold",family = "arial",colour = "black") ,axis.text.y = element_text(size = 15,face = "bold",family = "arial"),
                         plot.title = element_text(size=22, face="bold", hjust=0.5,vjust=0.4,family = "arial"))+
      ggtitle("Gene expression level")
  }

}
