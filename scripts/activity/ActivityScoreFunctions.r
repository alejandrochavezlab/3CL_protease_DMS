getMutations <- function() {
  c("*", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
    "S", "T", "V", "W", "Y")
}

makeMutNames <- function() {
  muts = getMutations()
  mutNames  = c()
  residues = 1:306
  for(resid in residues) {
    set = setForResidue[resid]
    refSite1 = refAtPos(resid-1, globalRef)
    refSite3 = refAtPos(resid+1, globalRef)
    
    mutNames = c(mutNames, paste0(set, "_residue", resid, "_", refSite1, muts, refSite3))
  }
  mutNames
}


ref_genetic <- function() { 
  seq_3CL = paste0("TACAAAATGTACAAAATGAGTGGTTTTAGAAAAATGGCATTCCCATCTGGTAAAGTTGAGGGTTGTATGGTACAAGTAACTTGTGGTACAACTACACTTAACGGTCTTTGGCTTGATGACGTAGTTTACTGTCCAAGACATGTGATCTGCACCTCTGAAGACATGCTTAACCCTAATTATGAAGATTTACTCATTCGTAAGTCTAATCATAATTTCTTGGTACAGGCTGGTAATGTTCAACTCAGGGTTATTGGACATTCTATGCAAAATTGTGTACTTAAGCTTAAGGTTGATACAGCCAATCCTAAGACACCTAAGTATAAGTTTGTTCGCATTCAACCAGGACAGACTTTTTCAGTGTTAGCTTGTTACAATGGTTCACCATCTGGTGTTTACCAATGTGCTATGAGGCCCAATTTCACTATTAAGGGTTCATTCCTTAATGGTTCATGTGGTAGTGTTGGTTTTAACATAGATTATGACTGTGTCTCTTTTTGTTACATGCACCATATGGAATTACCAACTGGAGTTCATGCTGGCACAGACTTAGAAGGTAACTTTTATGGACCTTTTGTTGACAGGCAAACAGCACAAGCAGCTGGTACGGACACAACTATTACAGTTAATGTTTTAGCTTGGTTGTACGCTGCTGTTATAAATGGAGACAGGTGGTTTCTCAATCGATTTACCACAACTCTTAATGACTTTAACCTTGTGGCTATGAAGTACAATTATGAACCTCTAACACAAGACCATGTTGACATACTAGGACCTCTTTCTGCTCAAACTGGAATTGCCGTTTTAGATATGTGTGCTTCATTAAAAGAATTACTGCAAAATGGTATGAATGGACGTACCATATTGGGTAGTGCTTTATTAGAAGATGAATTTACACCTTTTGATGTTGTTAGACAATGCTCAGGTGTTACTTTCCAATAA")
  seq_3CL
}

find_ref <- function(){
  dna = DNAStringSet(seq_3CL)
  ref= translate(dna, genetic.code=GENETIC_CODE, no.init.codon=TRUE)
  ref <- as.character(ref)
  ref <- substring(ref,6)
  return(ref)
}


filterRepeatedResidueFiles <- function(fnames) 
{
  res  = data.frame(t(sapply(fnames, parseFileName)))
  res$rep = as.numeric(res$rep)
  res$resid = as.numeric(res$resid)
  res$set = as.character(res$set)
  
  indexes = grep("R", res$set)
  setRs = res[indexes,]
  nonSetRs = res[-indexes,]
  
  rm_indexes = c()
  org_indexes = c()
  for(i in 1:nrow(setRs)) {
    resid = setRs[i, ]$resid
    rep = setRs[i, ]$rep
    set = setRs[i, ]$set
    
    index = which(nonSetRs$resid == resid & nonSetRs$rep == rep & nonSetRs$set != set)
    if(length(index) > 0) {
      rm_indexes = c(rm_indexes, index)
      org_indexes = c(org_indexes, i)
    }
  }
  

  nonSetRs = nonSetRs[-rm_indexes,]
  
  new_fnames = rownames(rbind(setRs, nonSetRs))
  new_fnames
}


makeFolderName <- function(base_folder, condition) {
  paste0(base_folder, condition, "/")
}

parseFileName <- function(fname) {
  tmp = strsplit(fname, "\\.")[[1]][1]
  splits = strsplit(tmp, "_")[[1]]
  
  list(set = splits[1], resid=substring(splits[2], 8), rep = substring(splits[3], 4))
}

readResidueFile <- function(folder, fname, condition, threshold, synCoding, remove_one_mismatch)
{
  full_fname = paste0(folder, fname)
  dat = read.csv(full_fname)
  
  if(remove_one_mismatch) {
    dat = dat[which(!dat$one_mismatch),] 
  }
  dat = dat[which(dat$count >= threshold), , drop=FALSE ]
  
  if(nrow(dat) <= 1) {
    return(NULL)
  }
  
  genetic_codes = apply(dat[, c("site_1", "site_2", "site_3")], 1, paste0, collapse='')
  names(genetic_codes) = NULL
  
  dna = DNAStringSet(genetic_codes)
  trans= as.character(translate(dna, genetic.code=GENETIC_CODE, no.init.codon=TRUE))
  
  parsedFname = parseFileName(fname)
  prefix = paste0(parsedFname$set, "_residue", parsedFname$resid, "_rep", parsedFname$rep)
  
  dat$id =   paste0(prefix, genetic_codes)
  dat$trans = trans
  
  mut = sapply(dat$trans, function(x){substr(x, 2, 2)})
  dat$mut = mut
  
  dat$set = parsedFname$set
  dat$condition = condition
  dat$rep = parsedFname$rep
  dat$gen_code = genetic_codes
  
  colnames(dat)[which(colnames(dat) == "resid")]  = "residue"
  
  native_wt_index = which(dat$WT)
  
  normFactor = dat[native_wt_index, "count"] + 1
  
  dat = dat[-native_wt_index,, drop=FALSE ]
  
  dat$count = (dat$count+1)/normFactor
  dat$corrected_count = (dat$corrected_count+1)/normFactor
  
  
  resid = dat$residue[1]
  dat$AA_WT = 0
  refSite2 = refAtPos(resid, globalRef)
  indexes = which(sapply(dat$trans, function(x){refSite2 == substr(x, 2, 2)}))
  if(length(indexes) == 0) {
  } else {
    dat[indexes, "AA_WT"] = 1
  }
  
  if(synCoding == FALSE) {
    dat = dat[, -which(colnames(dat) %in% c("id", "gen_code", "DistFromWT")) ]
    dat = as.data.frame(distinct(dat %>% group_by(mut) %>% 
                                   mutate(count=sum(count), corrected_count=sum(corrected_count))))
    dat$id = paste0(dat$set, "_residue", dat$residue, "_rep", dat$rep, "_", dat$trans)
  }
  
  dat
}


# whichRep: both, rep0, rep1
makeCountsDMS <- function(base_folder, condition, whichRep, threshold, synCoding, remove_one_mismatch) {
  folder = makeFolderName(base_folder, condition)
  fnames <-  list.files(folder, pattern="*.csv")
  if( whichRep %in% c("rep0", "rep1") ) {
    fnames = fnames[grep(whichRep,   fnames)]
  }
  
  fnames = filterRepeatedResidueFiles(fnames)
  
  ldf <- lapply(fnames, readResidueFile, folder = folder, condition=condition, 
                threshold=threshold, synCoding=synCoding, remove_one_mismatch=remove_one_mismatch)
  
  
  df = ldf[[1]]
  for(i in 2:length(ldf)) {
    df = rbind(df, ldf[[i]])  
  }
  df = df[order(as.numeric(df$residue)),]
  df
}


makeGluGal <- function(whichRep, gal_thr, glu_thr, normMethod, synCoding, remove_one_mismatch) {
  base_folder = "csv_files/"
  
  galDat = makeCountsDMS(base_folder, condition="Gal", whichRep = whichRep, threshold=gal_thr, synCoding=synCoding, remove_one_mismatch=remove_one_mismatch)
  gluDat = makeCountsDMS(base_folder, condition="Glu", whichRep = whichRep, threshold=glu_thr, synCoding=synCoding, remove_one_mismatch=remove_one_mismatch)
  
  
  sharedIDs = intersect(gluDat$id, galDat$id)
  
  glu = gluDat[match(sharedIDs, gluDat$id),]
  gal = galDat[match(sharedIDs, galDat$id),]
  all(glu$id == gal$id)
  
  glu_gal = gal
  
  glu_gal$condition = "Glu_Gal"
  if(normMethod == "ratio") {
    glu_gal$count =  log2(glu$count/gal$count)

  } else if(normMethod == "subtract") {
    glu_gal$count =  glu$count - gal$count
  }
  
  glu_gal
}

toWideFormat <- function(df) {
  sets = unique(df$set)
  df_all = list()
  for(set in sets) {
    df2 = df[which(df$set == set),]
    df2$id2 = paste0(df2$set, "_residue", df2$residue, "_", df2$trans)
    
    # counts in a wide format
    tmpWideDat = df2 %>% pivot_wider(names_from = id2, values_from = count)
    
    # remove non-count columns and transpose it
    indexes = intersect(grep("set", colnames(tmpWideDat)), grep("residue", colnames(tmpWideDat)))
    wideDat = t(tmpWideDat[, indexes])
    
    wideListDat = apply(wideDat, 1, function(x){ x[!is.na(x)]})
    # convert it to the list
    if(is.matrix(wideListDat)) {
      tmp = lapply(seq_len(ncol(wideListDat)), function(i) wideListDat[,i])
      names(tmp) = colnames(wideListDat)
      wideListDat = tmp
    }

    df_all = c(df_all, wideListDat)
  }
  
  mutNames = makeMutNames()
  res <- vector("list", length = length(mutNames))
  names(res) = mutNames
  
  indexes = match(names(df_all), mutNames)
  res[indexes] = df_all
  res
}

makeMutNames <- function() {
  muts = getMutations()
  mutNames  = c()
  residues = 1:306
  for(resid in residues) {
    set = setForResidue[resid]
    refSite1 = refAtPos(resid-1, globalRef)
    refSite3 = refAtPos(resid+1, globalRef)
    
    mutNames = c(mutNames, paste0(set, "_residue", resid, "_", refSite1, muts, refSite3))
  }
  mutNames
}


computeWT_counts <- function(df) {
  WT_all = df[which(df$AA_WT == 1), "count"]
  
  residues = 1:306
  
  WT_residue = sapply(residues, function(r) { 
    df[df$residue == r & df$AA_WT == 1, "count" ]
  } )
  names(WT_residue) = residues
  
  sets = sort(unique(df$set))
  WT_set = sapply(sets, function(s) { 
    df[df$set == s & df$AA_WT == 1, "count" ]
  } )
  names(WT_set) = sets
  
  list(WT_all=WT_all, WT_residue=WT_residue, WT_set=WT_set)  
}

annotateMuts <- function(mutIDs) {
  annots = data.frame()
  for(mut in mutIDs) {
    splits = strsplit(mut, "_")[[1]]
    
    annots = rbind(annots, data.frame(set = splits[1], resid=as.numeric(substring(splits[2], 8)), codon = splits[3], 
                                      mut = substr(splits[3], 2,2)))
  }
  
  annots$clinical_status = ""
  
  clinical_data <- read.csv("input_data/220325_3CLpro_clinical_variants.csv")
  for(i in 1:nrow(clinical_data)) {
    resid = clinical_data[i, ]$Residue
    mut = clinical_data[i, ]$mut_AA
    annots[which(annots$resid == resid &annots$mut == mut), "clinical_status"] = "clinical"
  }
  
  annots$WT = 0
  for(resid in 1:max(annots$resid)) {
    refSite = refAtPos(resid, globalRef)
    annots[which(annots$resid == resid & annots$mut == refSite), "WT"] = 1
  }
  
  annots
}

makeDMS_data <- function(whichRep, gal_thr, glu_thr, normMethod, synCoding, remove_one_mismatch){
  df = makeGluGal(whichRep, gal_thr = gal_thr, glu_thr = glu_thr, normMethod, synCoding, remove_one_mismatch)

  wide = toWideFormat(df)

  dmsWT = computeWT_counts(df)
  
  dms = list(counts=wide, wildType=dmsWT, annot=annotateMuts(names(wide)))
  class(dms) = "dms"
  dms
}

getSetsForResidue <- function() {
  base_folder = "csv_files/"
  condition = "Gal"
  folder = makeFolderName(base_folder, condition)
  fnames <-  list.files(folder, pattern="*.csv")
  fnames= fnames[grep( "rep0", fnames)]
  fnames = filterRepeatedResidueFiles(fnames)
  
  df = data.frame()
  for(fname in fnames) {
    df = rbind(df, unlist(parseFileName(fname) )  )
  }
  
  rownames(df)  = df$X.134.
  
  df2 = df[order(as.numeric(rownames(df))), 1, drop=FALSE]
  
  sets = df2[, 1]
  names(sets) = rownames(df2)
  setForResidue = sets
  setForResidue
}

#  access function for seq ref at site "pos"
refAtPos <- function(pos, ref) {
  substr(ref, pos + 1, pos + 1)
}

# WT_method: residue or set
computeAcitivityScores <- function(gal_thr, glu_thr, WT_method, whichRep, normMethod, synCoding, 
                                   remove_one_mismatch, onlyToWT) {
  dms_data = makeDMS_data(whichRep, gal_thr, glu_thr, normMethod, synCoding, remove_one_mismatch)
  
  result = data.frame()
  for(i in 1:length(dms_data$counts)) {
    mut_data = dms_data$counts[[i]]
    
    if(WT_method == "residue") {
      resid = dms_data$annot$resid[i]
      WT_data = dms_data$wildType$WT_residue[[resid]]
    } else if(WT_method == "set") {
      set = dms_data$annot$set[i]
      WT_data = dms_data$wildType$WT_set[[set]]
    } else {
      stop("unkown WT_method")
    }
    if(length(mut_data) > 1 & length(WT_data) > 1) {
      ttest_res = t.test(mut_data, WT_data, alternative="less")
      result = rbind(result, c(ttest_res$statistic, ttest_res$p.value, length(mut_data), length(WT_data)))
    } else {
      result = rbind(result, c(NA, NA, NA, NA))
    }
  }
  colnames(result) = c("AS", "AS_pvalue", "nr_mut", "nr_wt")
  result$AS_fdr = p.adjust(result$AS_pvalue, method = "fdr")
  rownames(result) = names(dms_data$counts)
  result = cbind(result, dms_data$annot)
  
  result = result[, c("nr_mut", "nr_wt", "AS", "AS_pvalue", "AS_fdr", "set",  "resid", "codon", 
                      "mut",  "WT", "clinical_status")]
  result
}

.rescaleToWT_STOP <- function(values, stop_value, wt_value){
  (values - stop_value) /(wt_value - stop_value) -1
}

rescaleActivityScores <- function(act, onlyToWT) {
  act$AS_raw = act$AS
  
  sets = unique(act$set)
  for(s in sets) {
    indexes = which(act$set == s)
    t = act[indexes, ]
    
    stop_value = mean(t$AS[t$mut == "*"],  na.rm=T )
    wt_value = mean(t$AS[t$WT == 1],  na.rm=T )
    
    if(onlyToWT) {
      a = t$AS - wt_value
    } else {
      a = .rescaleToWT_STOP(t$AS, stop_value, wt_value)  
    }

    act[indexes, "AS_scaled"] = a
  }
  
  act$AS = act$AS_scaled
  act = act[, c("nr_mut", "nr_wt", "AS_raw", "AS", "AS_pvalue", "AS_fdr", "set",  "resid", "codon",
                "mut",  "WT", "clinical_status")]
  act
}
