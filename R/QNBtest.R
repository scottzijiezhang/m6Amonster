#' @title QNBtest
#' @param x The monster data object
#' @export
QNBtest <- function(x){
  library(QNB)
  Ctl_id <- which( x$X %in%  unique(x$X)[1] )
  Treat_id <-which( x$X %in%  unique(x$X)[2] )
  QNB.res <- QNB::qnbtest(control_ip = x$jointPeak_ip[,Ctl_id],
                            treated_ip = x$jointPeak_ip[,Treat_id],
                            control_input = x$jointPeak_input[,Ctl_id],
                            treated_input = x$jointPeak_input[,Treat_id],plot.dispersion = FALSE)
  colnames(QNB.res) <- c("p.treated","p.control", "log2.RR" ,"beta",  "p_value", "q","padj"    )

  result <- c(x,list("all.est.peak" = QNB.res))
  return(result)
}
