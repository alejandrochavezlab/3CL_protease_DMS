rm(list=ls())
library(foreach)
library('stringr')
library(readr)
library(ggplot2)
library(ggExtra)
library(Biostrings)

source("scripts/activity/ActivityScoreFunctions.r")


setForResidue = getSetsForResidue()

seq_3CL = ref_genetic()
globalRef = find_ref()

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

Gal_thr = Glu_thr = c(0, 10, 30, 50, 100)
z = c(FALSE, TRUE)
settings = expand.grid(Glu_thr=Glu_thr, Gal_thr=Gal_thr, remove_one_mismatch=c(FALSE, TRUE), 
                         synCoding=c(TRUE), WT_method = c("set"))

t0 = proc.time()
corr_results = foreach(i=1:nrow(settings), .combine='rbind') %dopar% {
  library(Biostrings)
  library(data.table)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  
  
  glu_thr = settings[i, "Glu_thr"]
  gal_thr = settings[i, "Gal_thr"]
  remove_one_mismatch = settings[i, "remove_one_mismatch"]
  synCoding = settings[i, "synCoding"]
  WT_method = as.character(settings[i, "WT_method"])
  
  
  act = computeAcitivityScores(gal_thr = gal_thr, glu_thr = glu_thr, WT_method = WT_method, 
                               whichRep = "both", normMethod = "ratio", synCoding=synCoding, 
                               remove_one_mismatch = remove_one_mismatch)
  missing = mean(is.na(act$AS_pvalue))
  # 
  act0 = computeAcitivityScores(gal_thr = gal_thr, glu_thr = glu_thr, WT_method = WT_method, 
                                whichRep = "rep0", normMethod = "ratio", synCoding=synCoding, 
                                remove_one_mismatch = remove_one_mismatch)
  
  act1 = computeAcitivityScores(gal_thr = gal_thr, glu_thr = glu_thr, WT_method = WT_method, 
                                whichRep = "rep1", normMethod = "ratio", synCoding=synCoding, 
                                remove_one_mismatch = remove_one_mismatch)
  
  corr = cor(act0$AS, act1$AS, method = "spearman", use="pairwise.complete.obs")
  
  print(c(corr, missing))
  print(proc.time() - t0)
  c(corr, missing)
}

print(proc.time() - t0)
colnames(corr_results) = c("corr", "missing")
corr_results = cbind(settings[, 1:3], corr_results[, c("corr", "missing")])
rownames(corr_results) = NULL
write.csv(corr_results, file="outputs/dms_correlation.csv", row.names = FALSE)

