rm(list=ls())
source("scripts/activity/ActivityScoreFunctions.r")
library(tidyverse)
library('Biostrings')
library(ggExtra)
source("scripts/activity/plotting.r")


runActivityScript <- function(gal_thr, glu_thr, WT_method = "set", normMethod = "ratio", remove_one_mismatch = TRUE,
                              synCoding = TRUE, output_folder, onlyToWT){

  act = computeAcitivityScores(gal_thr = gal_thr, glu_thr = glu_thr, WT_method = WT_method, 
                               whichRep = "both", normMethod = normMethod, synCoding=synCoding,
                               remove_one_mismatch = remove_one_mismatch)
  
  act0 = computeAcitivityScores(gal_thr = gal_thr, glu_thr = glu_thr, WT_method = WT_method,
                                whichRep = "rep0", normMethod = normMethod, synCoding=synCoding,
                                remove_one_mismatch = remove_one_mismatch)
  
  act1 = computeAcitivityScores(gal_thr = gal_thr, glu_thr = glu_thr, WT_method = WT_method,
                                whichRep = "rep1", normMethod = normMethod, synCoding==synCoding,
                                remove_one_mismatch = remove_one_mismatch)
  
  indexes = which(is.na(act0$AS) & is.na(act1$AS))
  act[indexes,] = act0[indexes,] = act1[indexes,] = NA
  
  act = rescaleActivityScores(act, onlyToWT)
  act0 = rescaleActivityScores(act0, onlyToWT)
  act1 = rescaleActivityScores(act1, onlyToWT)
  
    # Task 1: plot histogram for syn coding
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  
  pdf(paste0(output_folder, "histogram.pdf"))
  hist.plot(act)
  dev.off()
  
  write.csv(act, file=paste0(output_folder, "activity_scores.csv"))
  
  # Task 2: correlation/scatter plot
  
  spearman_cor = cor(act0$AS, act1$AS, method = "spearman", use="pairwise.complete.obs")
  pearson_cor = cor(act0$AS, act1$AS, method = "pearson", use="pairwise.complete.obs")
  sum(is.na(act$AS))
  
  pdf(paste0(output_folder, "scatter_plot.pdf"))
  print(scattor.plot(act0, act1))
  dev.off()
  
  write.csv(act0, file=paste0(output_folder, "activity_scores_rep0.csv"))
  write.csv(act1, file=paste0(output_folder, "activity_scores_rep1.csv"))
  
  act0$AS[which(act0$AS< -10)]  = -10
  act1$AS[which(act1$AS< -10)]  = -10
  
  
  pdf(paste0(output_folder, "scatter_truncated.pdf"))
  print(scattor.plot(act0, act1))
  dev.off()
  c(spearman_cor, pearson_cor)
}

setForResidue = getSetsForResidue()
seq_3CL = ref_genetic()
globalRef = find_ref()


# set the final settings
gal_thr = 0
glu_thr = 11  # equivalent to gal_thr >= 10
output_folder = "outputs/results/normalized_to_wt_and_stop/"
runActivityScript(gal_thr, glu_thr, output_folder=output_folder, onlyToWT=FALSE)

# output_folder = "outputs/results/normalized_only_to_wt/"
# runActivityScript(gal_thr, glu_thr, output_folder=output_folder, onlyToWT=TRUE)
