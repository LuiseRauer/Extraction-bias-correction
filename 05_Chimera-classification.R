################################################################################
#
# Zymo-Qiagen comparison mock analysis
# Chimera identification by longest common substring (LCS)
#
################################################################################

# Load required packages
library(tidyverse)
library(PTXQC) # LCS function

################################################################################
#
# Longest common substring for chimera detection
#
################################################################################

# ------------------------------------------------------------------------------
# Chimera identification algorithm
# ------------------------------------------------------------------------------

# Detect chimeras
Chim_detect <- function(Seq, Comparison_seqs) {
  # LCS round 1
  Res1 <- sapply(Comparison_seqs, function(x) LCS(Seq, x))
  LCS_max1 <- which.max(nchar(Res1))
  Seq2 <- sub(Res1[LCS_max1], "", Seq)
  # LCS round 2
  if(nchar(Seq2) > 7) {
    Res2 <- sapply(Comparison_seqs, function(x) LCS(Seq2, x))
    LCS_max2 <- which.max(nchar(Res2))
    Seq3 <- sub(Res2[LCS_max2], "", Seq2)
    # LCS round 3
    if(nchar(Seq3) > 7) {
      Res3 <- sapply(Comparison_seqs, function(x) LCS(Seq3, x))
      LCS_max3 <- which.max(nchar(Res3))
      Seq4 <- sub(Res3[LCS_max3], "", Seq3)
      # LCS round 4
      if(nchar(Seq4) > 7) {
        Res4 <- sapply(Comparison_seqs, function(x) LCS(Seq4, x))
        LCS_max4 <- which.max(nchar(Res4))
      }
    }
  }
  return(list(
    # Original sequence:
    Seq, 
    # Remaining sequence after removing the LCS
    Seq2, 
    # Same for remaining sequences after removing more LCS
    if(nchar(Seq2) > 7) Seq3 else '',
    if(nchar(Seq2) > 7 && nchar(Seq3) > 7) Seq4 else '',
    # Number in Table_16S_279_unique with the highest LCS
    LCS_max1, 
    # Remaining numbers of the highest LCS (not needed later)
    if(nchar(Seq2) > 7) LCS_max2 else '', 
    if(nchar(Seq2) > 7 && nchar(Seq3) > 7) LCS_max3 else '',
    if(nchar(Seq2) > 7 && nchar(Seq3) > 7 && nchar(Seq4) > 7) LCS_max4 else '',
    # Complete export of the LCS results per step (length each: 21)
    res1, 
    if(nchar(Seq2) > 7) res2, 
    if(nchar(Seq2) > 7 && nchar(Seq3) > 7) res3,
    if(nchar(Seq2) > 7 && nchar(Seq3) > 7 && nchar(Seq4) > 7) res4))
}

Sys.time() # Takes 26s for 10 Sequences (~ 2h30min for 3391 ASVs)
CCC_all <- lapply(
  seqtab.final[, "Sequence"], function(x)
    Chim_detect(x, Table_16S_279_unique$Sequence))
Sys.time()
###save(CCC_all, file = "CCC_all_2024-08-20.RData")
###load("CCC_all_2024-08-20.RData", verbose = TRUE)

# ------------------------------------------------------------------------------
# Extract the "annotated" (remaining) sequence per LCS step
# ------------------------------------------------------------------------------

# Extract list items 1-4 per ASV
CCC_seq <- lapply(CCC_all, function(x) unlist(x[c(1:4)]))
# Transform list to dataframe
CCC_seq <- data.frame(t(do.call(cbind, CCC_seq)))
# Replace empty sequences by NA
CCC_seq[CCC_seq == ""] <- NA
# Set colnames
colnames(CCC_seq) <- c("Seq_L1", "Seq_L2", "Seq_L3", "Seq_L4")

# ------------------------------------------------------------------------------
# Extract annotation/species name per LCS step
# ------------------------------------------------------------------------------

# Extract 'Res1', the complete LCS result of the first round
CCC_LCS_L1 <- lapply(CCC_all, function(x) 
  unlist(lapply(x[9], function(y) nchar(y))))
# Extract the top annotation number per sequence (keeps ties) 
CCC_anno_L1 <- lapply(CCC_LCS_L1, function(x) 
  which(x == unique(sort(x, decreasing = TRUE))[1]))
# Map the sequence number with the actual reference species name (keeps ties)
CCC_anno_L1 <- lapply(CCC_anno_L1, function(x) 
  unique(Table_16S_279_unique[x, "ASV.ID"]))

# Check that ambiguous annotations only occur with Esch./Salm. or low LCS values
intersect(intersect(
  which(lapply(CCC_anno_L1, length) > 1), # Sequence numbers with >1 annotation
  which(unlist(lapply(CCC_LCS_L1, max)) > 50)), # Sequence numbers with small LCS
  which(lapply(CCC_anno_L1, function(x) all(grepl("Esch|Salm", x))) == FALSE))
# No sequences left, so I can take the first species only

# Extract only the first annotation
CCC_anno_L1 <- data.frame(Taxon_L1 = sapply(CCC_anno_L1, "[[", 1))

# Extract 'Res2', the complete LCS result of the second round
CCC_LCS_L2 <- lapply(CCC_all, function(x) 
  unlist(lapply(x[10], function(y) nchar(y))))
# Extract the top annotation number per sequence (keeps ties) 
CCC_anno_L2 <- lapply(CCC_LCS_L2, function(x) 
  which(x == unique(sort(x, decreasing = TRUE))[1]))
# Map the sequence number with the actual reference species name (keeps ties)
CCC_anno_L2 <- lapply(CCC_anno_L2, function(x) 
  unique(Table_16S_279_unique[x, "ASV.ID"]))

# Check that ambiguous annotations only occur with Esch./Salm. or low LCS values
intersect(intersect(
  which(lapply(CCC_anno_L2, length) > 1), # Sequence numbers with >1 annotation
  which(unlist(lapply(CCC_LCS_L2, max)) > 50)), # Sequence numbers with small LCS
  which(lapply(CCC_anno_L2, function(x) all(grepl("Esch|Salm", x))) == FALSE))
# No sequences left, so I can take the first species only

# Extract only the first annotation
CCC_anno_L2 <- lapply(CCC_anno_L2, function(x) replace(x, length(x) == 0, NA))
CCC_anno_L2 <- data.frame(Taxon_L2 = sapply(CCC_anno_L2, "[[", 1))

# Extract 'Res3', the complete LCS result of the second round
CCC_LCS_L3 <- lapply(CCC_all, function(x) 
  unlist(lapply(x[11], function(y) nchar(y))))
# Map the sequence number with the actual reference species name (keeps ties)
CCC_anno_L3 <- lapply(CCC_LCS_L3, function(x) 
  which(x == unique(sort(x, decreasing = TRUE))[1]))
# Map the sequence number with the actual reference species name (keeps ties)
CCC_anno_L3 <- lapply(CCC_anno_L3, function(x) 
  unique(Table_16S_279_unique[x, "ASV.ID"]))

# Check that ambiguous annotations only occur with Esch./Salm. or low LCS values
intersect(intersect(
  which(lapply(CCC_anno_L3, length) > 1), # Sequence numbers with >1 annotation
  which(unlist(lapply(CCC_LCS_L3, max)) > 50)), # Sequence numbers with small LCS
  which(lapply(CCC_anno_L3, function(x) all(grepl("Esch|Salm", x))) == FALSE))
# No sequences left, so I can take the first species only

# Extract only the first annotation
CCC_anno_L3 <- lapply(CCC_anno_L3, function(x) replace(x, length(x) == 0, NA))
CCC_anno_L3 <- data.frame(Taxon_L3 = sapply(CCC_anno_L3, "[[", 1))

# Extract 'Res4', the complete LCS result of the second round
CCC_LCS_L4 <- lapply(CCC_all, function(x) 
  unlist(lapply(x[12], function(y) nchar(y))))
# Map the sequence number with the actual reference species name (keeps ties)
CCC_anno_L4 <- lapply(CCC_LCS_L4, function(x) 
  which(x == unique(sort(x, decreasing = TRUE))[1]))
# Map the sequence number with the actual reference species name (keeps ties)
CCC_anno_L4 <- lapply(CCC_anno_L4, function(x) 
  unique(Table_16S_279_unique[x, "ASV.ID"]))

# Check that ambiguous annotations only occur with Esch./Salm. or low LCS values
intersect(intersect(
  which(lapply(CCC_anno_L4, length) > 1), # Sequence numbers with >1 annotation
  which(unlist(lapply(CCC_LCS_L4, max)) > 50)), # Sequence numbers with small LCS
  which(lapply(CCC_anno_L4, function(x) all(grepl("Esch|Salm", x))) == FALSE))
# No sequences left, so I can take the first species only

# Extract only the first annotation
CCC_anno_L4 <- lapply(CCC_anno_L4, function(x) replace(x, length(x) == 0, NA))
CCC_anno_L4 <- data.frame(Taxon_L4 = sapply(CCC_anno_L4, "[[", 1))

# ------------------------------------------------------------------------------
# Extract length of the highest LCS match
# ------------------------------------------------------------------------------

# Extract the top three annotation values of round L1
CCC_len_L1 <- data.frame(Len_L1 = sapply(CCC_LCS_L1, function(x) 
    unique(sort(x, decreasing = TRUE))[1]))

# ------------------------------------------------------------------------------
# Combine all data
# ------------------------------------------------------------------------------

# Make a table out of that...
CCC_res <- do.call("cbind.data.frame", list(
  CCC_seq, nchar = apply(CCC_seq, 2, nchar),
  CCC_len_L1, CCC_anno_L1, CCC_anno_L2, CCC_anno_L3, CCC_anno_L4))
str(CCC_res)

# Add ASV table and RDP/LV annotation
CCC_res <- merge(
  CCC_res, seqtab.final, by.x = "Seq_L1", by.y = "Sequence", all = TRUE)
# Add minimum LV distance + annotation
CCC_res <- merge(
  CCC_res, Min_sequ_dist, by.x = "ASV_ID", by.y = "ID", all = TRUE)
# Add top 2/3 annotation if difference in matched length is < 10?

################################################################################
#
# Classification into sequence errors, chimeras and Unknowns
#
################################################################################

# 97% identity: sequence error, 95% identity with two or three parts: chimera
CCC_res <- CCC_res %>% 
  mutate(Category = case_when(
    Len_L1 == 279 ~ "Exact_match",
    Min_LV_dist < 9 ~ "Amplification_error",
    (nchar.Seq_L3 <= 13 | is.na(nchar.Seq_L3) | nchar.Seq_L4 <= 13 |
       is.na(nchar.Seq_L4)) ~ "Mock_chimera",
    TRUE ~ "Unknown"))

CCC_list <- list(
  CCC_anno_L1, CCC_anno_L2, CCC_anno_L3, CCC_anno_L4, CCC_LCS_L1, CCC_LCS_L2,
  CCC_LCS_L3, CCC_LCS_L4, CCC_len_L1, CCC_seq, CCC_res)
###save(CCC_list, file = "R_QZmock_CCC_res_parts_2023-08-21.RData")
###load("R_QZmock_CCC_res_parts_2023-08-21.RData", verbose = TRUE)

rm(CCC_anno_L1, CCC_anno_L2, CCC_anno_L3, CCC_anno_L4, CCC_LCS_L1, CCC_LCS_L2,
   CCC_LCS_L3, CCC_LCS_L4, CCC_len_L1, CCC_list, CCC_seq, Chim_detect)  

# R objects after running this script:
# "CCC_all", "CCC_res"
