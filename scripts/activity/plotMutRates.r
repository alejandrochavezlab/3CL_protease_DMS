rm(list=ls())
library(Biostrings)
library(dplyr)
source("scripts/activity/ActivityScoreFunctions.r")

rate.plot <- function(mut_data){
  data = t(as.matrix(mut_data[,c(2,3,4)]))
  colnames(data) <- seq(1,306)
  data = apply(data, 2, function(x){x/sum(x)})
  barplot(data, 
          col=c("#4085ed","#cc195d","#f7c200"),
          border="white", 
          font.axis=2,las = 1, 
          xlab="Residue", ylab = "Fraction")
  legend("bottom", inset=c(-0.2, 0),legend=c("Non-synonymous","Synonymous","Stop"), fill =c("#4085ed","#cc195d","#f7c200"), cex=0.8,
         box.lty=0, xpd = TRUE,  bty = "n")
}


computeMutationRatesPerResidue <- function(dmsData) {

  dmsData$mut2 = dmsData$mut
  dmsData[which(dmsData$mut == "*"), "mut2"] = "stop"
  dmsData[which(dmsData$mut != "*"), "mut2"] = "Non_syn"
  dmsData[which(dmsData$AA_WT==1), "mut2"] = "Syn"
  
  table(dmsData$mut2)
  
  
  df2 = dmsData %>% group_by(residue, mut2, rep) %>%  summarise(count=sum(count))
  df2$mut2 = as.factor(df2$mut2)
  df2 = as.data.frame(df2)
  
  a = reshape(df2, timevar = "mut2", idvar = c("residue", "rep"), direction = "wide")
  a[is.na(a)] = 0
  
  a[, 3:5] = a[, 3:5] / apply(a[, 3:5], 1, sum)
  rownames(a) = NULL
  
  b = a %>% group_by(residue) %>%  summarise(count.Non_syn=mean(count.Non_syn), count.stop = mean(count.stop), 
                                             count.Syn=mean(count.Syn) )
  
  b = as.data.frame(b)
  b
}

figure_folder = "outputs/results/mutation_rates/"
dir.create(figure_folder, recursive = TRUE, showWarnings = FALSE)

setForResidue = getSetsForResidue()

seq_3CL = ref_genetic()
globalRef = find_ref()

base_folder = "csv_files/"
whichRep = "both"
remove_one_mismatch = TRUE
synCoding = TRUE



gal_thr = 0
galDat = makeCountsDMS(base_folder, condition="Gal", whichRep = whichRep, threshold=gal_thr, synCoding=TRUE, remove_one_mismatch=remove_one_mismatch)
mutRates = computeMutationRatesPerResidue(galDat)

colnames(mutRates)[2:4] = c("MR_Non_syn", "MR_stop", "MR_Syn")

write.csv(mutRates, file=paste0(figure_folder, "mut_rates_gal.csv"), row.names = F)
pdf(paste0(figure_folder, "mut_rates_gal.pdf"), width = 9, height = 6)
rate.plot(mutRates)
dev.off()

glu_thr = 11
gluDat = makeCountsDMS(base_folder, condition="Glu", whichRep = whichRep, threshold=glu_thr, synCoding=TRUE, remove_one_mismatch=remove_one_mismatch)
mutRates = computeMutationRatesPerResidue(gluDat)

colnames(mutRates)[2:4] = c("MR_Non_syn", "MR_stop", "MR_Syn")

write.csv(mutRates, file=paste0(figure_folder, "mut_rates_glu_10.csv"), row.names = F)
pdf(paste0(figure_folder, "mut_rates_glu_10.pdf"), width = 9, height = 6)
rate.plot(mutRates)
dev.off()

