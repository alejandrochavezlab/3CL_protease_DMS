

scattor.plot <- function(data_act0,data_act1){
  data = as.data.frame(cbind("AS0"=data_act0$AS,"AS1"=data_act1$AS))
  r = round(cor(data$AS0, data$AS1, method = "spearman", use="pairwise.complete.obs"),digits=2)
  
  p <- ggplot(data, aes(x=AS0, y=AS1)) +
    geom_point(color="#00AFBB") +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    labs(y="Activity score in rep1", x = "Activity score in rep0")
  p <- p + ggplot2::annotate("text", x = -7, y = 7, label = paste0("italic(Correlation) == ",r),parse = TRUE)
  
  # with marginal histogram
  ggMarginal(p, type="histogram",xparams = list(fill = "#E7B800",na.rm=TRUE),yparams = list(fill = "#FC4E07",na.rm=TRUE)) 
}

hist.plot <- function(act_data, scaled=FALSE){
  AS_var = "AS"  
  if( scaled) {
    AS_var = "AS_scaled"
  }
  all_res = act_data[, AS_var]
  clin_res = (act_data[which(act_data$clinical_status =="clinical"), AS_var])
  stop_res = (act_data[which(act_data$mut == "*"), AS_var])
  #plot abline:
  nrResidues = length(unique(rownames(act_data)))
  alpha =  0.05
  y_max = 2
  
  avgSampleSize = round(mean(act_data$nr_mut + act_data$nr_wt, na.rm=TRUE))
  R = qt(alpha/nrResidues, avgSampleSize)
  hist(all_res, breaks=50, xlim=c(-3,3),ylim=c(0,y_max), col=rgb(1,0,0,0.5), xlab="Activity score", ylab="Frequency", freq = FALSE, main = "")
  hist(clin_res, breaks=30, xlim=c(-3,3), col=rgb(0,0,1,0.5),ylim=c(0,y_max), add=T,freq = FALSE)
  hist(stop_res, breaks=30, xlim=c(-3,3),ylim=c(0,y_max), col=rgb(0.2,0.8,0.5,0.5), add=T,freq = FALSE)
  legend("topright", inset=c(0,0), legend=c("All residues","Clinical isolates","Stop residues"), col=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5),rgb(0.2,0.8,0.5,0.5)), pt.cex=1, pch=15 )
  
  abline(v=R, col=c("#f3dc95"), lty=c(2), lwd=c(3))
}
