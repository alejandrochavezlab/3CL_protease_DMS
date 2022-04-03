rm(list=ls())
library(EnvStats)

act0 = read.csv("outputs/results/normalized_to_wt_and_stop/activity_scores_rep0.csv")
act1 = read.csv("outputs/results/normalized_to_wt_and_stop/activity_scores_rep1.csv")

nr_muts = apply(cbind(act0$nr_mut, act1$nr_mut), 1, max)

x = nr_muts[which(!is.na(nr_muts))] 

length(x)

output_folder = "outputs/results/number_of_codings/"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)


pdf(paste0(output_folder,"number_of_codings.pdf"))
ecdfPlot(x, discrete=TRUE, main ='', ylab="Fraction", xlab="Number of codings",
         ecdf.col = "blue")
dev.off()

a = ecdfPlot(x, discrete=TRUE, main ='', ylab="Fraction", xlab="Number of codings",
         ecdf.col = "blue")

dat = cbind(a$Order.Statistics, a$Cumulative.Probabilities)
colnames(dat) = c("Number of codings", "Cumulative Probability")

write.csv(dat, file=paste0(output_folder,"number_of_codings.csv"))
