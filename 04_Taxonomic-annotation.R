################################################################################
#
# Zymo-Qiagen comparison mock analysis
# Taxonomic annotation by database and exact matching with exp. sequences
#
################################################################################

# Load required packages
library(dada2) # assignTaxonomy
library(stringdist) # stringdistmatrix
library(ape) # dist.dna
library(tidyverse) # case_when in Exact_grepl_function
library(reshape2) # melt
library(ComplexHeatmap) # Heatmap of expected sequence distance

################################################################################
#
# Taxonomic annotation by databases: RDP
#
################################################################################

# Downloaded the file RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz from 
#   https://zenodo.org/record/4735821#.YO9BsUHgo2w and saved in output
#   directory (downloaded 2021-07-14)

print(Sys.time()) # takes approx. 1h10min
set.seed(1)
taxa_RDP <- assignTaxonomy(
  colnames(seqtab), 
  paste0(path_output, "/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz"),
  multithread = FALSE)
print(Sys.time())
###save(taxa_RDP, file = "RDP_taxa_forward_2023-05-24.RData")
###load("RDP_taxa_forward_2023-05-24.RData", verbose = TRUE)

# Add RDP annotation to seqtab.final
seqtab.final <- merge(seqtab.final, taxa_RDP, by.x = "Sequence", by.y = 0)
rm(taxa_RDP)

################################################################################
#
# Load and check expected 16S sequences for exact matching
#
################################################################################

# ------------------------------------------------------------------------------
# Load 16S expected WGS data for mock samples
# ------------------------------------------------------------------------------

# List all 16S files
Files_16S_mock <- list.files(
  paste0(path_proj, "Input/ZymoBIOMICS.STD.refseq.v2/ZymoBIOMICS.STD.refseq.v2/ssrRNAs"), 
  pattern = ".fasta", full.names = TRUE, recursive = TRUE)
# Load 16S files and store in a list
List_16S_mock <- lapply(Files_16S_mock, function(x) read.table(x))

# Reformat FASTA files to a table
Table_16S_mock <- lapply(List_16S_mock, function(x) 
  data.frame(ASV.ID = substring(x$V1[seq(1, length(x$V1), 2)], first = 2),
             Sequence = x$V1[seq(1, length(x$V1), 2)+1]))
Table_16S_mock <- Reduce(rbind.data.frame, Table_16S_mock)
rm(Files_16S_mock, List_16S_mock)

# Remove fungi
Table_16S_mock <- Table_16S_mock[!grepl("Crypt|Sacch", Table_16S_mock$ASV.ID), ]

# ------------------------------------------------------------------------------
# Load 16S expected WGS data for spike-in samples
# ------------------------------------------------------------------------------

# List all 16S files
Files_16S_spike <- list.files(
  paste0(path_proj, "Input/D6321.refseq/16S"), 
  pattern = ".fasta", full.names = TRUE, recursive = TRUE)
# Load 16S files in a loop and store in a list
List_16S_spike <- lapply(Files_16S_spike, function(x) read.table(x))

# Reformat FASTA files to a table
Table_16S_spike <- lapply(List_16S_spike, function(x) 
  data.frame(ASV.ID = substring(x$V1[seq(1, length(x$V1), 2)], first = 2),
             Sequence = x$V1[seq(1, length(x$V1), 2)+1]))
Table_16S_spike <- Reduce(rbind.data.frame, Table_16S_spike)
rm(Files_16S_spike, List_16S_spike)

# Remove consensus row
Table_16S_spike <- Table_16S_spike[Table_16S_spike$ASV.ID != "Consensus", ]

# ------------------------------------------------------------------------------
# Combine and check expected sequences
# ------------------------------------------------------------------------------

# Rbind mock and spike expected sequences
Table_16S <- rbind.data.frame(Table_16S_mock, Table_16S_spike)
rm(Table_16S_mock, Table_16S_spike)

# Remove copy information from species/gene name (because of counting errors)
Table_16S$ASV.ID <- sub("([A-Za-z]+_[A-Za-z]+).*", "\\1", Table_16S$ASV.ID)

# Add a new correct numbering of 16S copies
Table_16S <- Table_16S %>%
  group_by(ASV.ID) %>%
  mutate(ASV.ID_long = paste0(ASV.ID, "_", row_number())) %>% 
  as.data.frame()

# Check for ambiguous bases
summary(sapply(strsplit(Table_16S$Sequence, split = ""), function(x)
  which(!x %in% c("A", "C", "G", "T"))))
# No ambiguous bases left (only in the two fungi)

################################################################################
#
# Exact taxonomy with tolerance (Levenshtein distance)
#
################################################################################

# ------------------------------------------------------------------------------
# Extract forward subregion of the V1-V3 HVR from 16S genes
# ------------------------------------------------------------------------------

# Check if the primer can be detected in every expected sequence (4x FALSE)
table(grepl("AGAGTTTGAT(C|T)(A|C)TGGCTCAG", Table_16S$Sequence))
# For which taxa does the primer not match? 
Table_16S[!grepl("AGAGTTTGAT(C|T)(A|C)TGGCTCAG", Table_16S$Sequence), 1]
# Primers match, except for four P. aeruginosa sequences.
# Check if the four P. aeruginosa sequences are identical - yes.
stringdistmatrix(setNames(
  Table_16S[Table_16S$ASV.ID == "Pseudomonas_aeruginosa", "Sequence"], 
  Table_16S[Table_16S$ASV.ID == "Pseudomonas_aeruginosa", "ASV.ID_long"]), 
  method = "lv", useNames = "names")
# The primer in P. aeruginosa misses the first two bases.

# But the "shorter" primer is present in all 60 sequences
table(grepl("AGTTTGAT(C|T)(A|C)TGGCTCAG", Table_16S$Sequence))
# The shorter primer is only detected max. once per expected sequence 
stringr::str_count(string = Table_16S$Sequence, 
                   pattern = "AGTTTGAT(C|T)(A|C)TGGCTCAG")

# Check if the rev. complement of the reverse primer is present in all bacteria
table(grepl(rev_comp("ATTACCGCGGCTGCTGG"), Table_16S$Sequence))
# Present in all bacteria

# Cut sequences after the shortened forward primer
Table_16S_279 <- Table_16S
Table_16S_279$Sequence <- gsub(".*AGTTTGAT(C|T)(A|C)TGGCTCAG", "", 
                               Table_16S_279$Sequence)
# Cut sequences after 279 bases (reverse primer not covered with forward reads)
Table_16S_279$Sequence <- sapply(Table_16S_279$Sequence, function(x) 
  paste0(strsplit(x, split = "")[[1]][1:279], collapse = ""))
# Check length of expected sequences - 279 
table(nchar(Table_16S_279$Sequence))

# ------------------------------------------------------------------------------
# Levenshtein distance between expected sequences
# ------------------------------------------------------------------------------

# Calculate expected Levenshtein distance between cut 16S sequences
Exp_dist_279 <- as.matrix(stringdistmatrix(
  setNames(Table_16S_279$Sequence, Table_16S_279$ASV.ID),
  method = "lv", useNames = "names"))
# Empty the lower half of the distance matrix to not double values when melting
Exp_dist_279[lower.tri(Exp_dist_279)] <- NA

# Summarise mean LV distance between species
Exp_dist_279_mean <- Exp_dist_279 %>% 
  # Create long format
  melt(., na.rm = TRUE, varnames = c("Species1", "Species2"), value.name = "LV") %>% 
  # Calculate mean distance per species
  group_by(Species1, Species2) %>% summarise(LV = mean(LV)) %>%
  # Create wide format
  pivot_wider(names_from = Species2, values_from = LV) %>% 
  column_to_rownames("Species1") %>%
  # Save as matrix for next step
  as.matrix()
# Fill lower triangle again with mean values for visualisation
Exp_dist_279_mean[lower.tri(Exp_dist_279_mean)] <- 
  t(Exp_dist_279_mean)[lower.tri(Exp_dist_279_mean)]
# Replace underscore by space
row.names(Exp_dist_279_mean) <- gsub("_", " ", row.names(Exp_dist_279_mean))
colnames(Exp_dist_279_mean) <- gsub("_", " ", colnames(Exp_dist_279_mean))

# Plot the heatmap
svg("Output/Plots/S03_Expected-LV-distance_2024-11-04.svg", 
    width = 4.8, height = 4.8)
htmp <- Heatmap(
  Exp_dist_279_mean, name = "Levenshtein distance",
  # Write value per cell
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf(
      if(Exp_dist_279_mean[i, j] > 10 | Exp_dist_279_mean[i, j] == 0) "%.0f" 
      else "%.1f", Exp_dist_279_mean[i, j]), x, y, gp = gpar(fontsize = 8))},
  col = c("white", RColorBrewer::brewer.pal(name = "YlGnBu", n = 9)[1:7]),
  column_names_gp = gpar(fontsize = 9, fontface = "italic"), 
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  show_heatmap_legend = FALSE)
# Plot the heatmap
htmp = draw(htmp)
# Extract the (hidden) legend from the plot
lgd = color_mapping_legend(htmp@ht_list[[1]]@matrix_color_mapping, plot = FALSE)
# Draw the legend in the plot
draw(lgd, x = unit(10.1, "cm"), y = unit(2.3, "cm"))
dev.off()
rm(htmp, lgd)

# ------------------------------------------------------------------------------
# Levenshtein distance between observed and expected sequences
# ------------------------------------------------------------------------------

# Only keep unique sequences
Table_16S_279_unique <- 
  Table_16S_279[!duplicated(Table_16S_279[, c("ASV.ID", "Sequence")]), ]
row.names(Table_16S_279_unique) <- NULL
# Levenshtein distance between expected cut 16S and observed forward sequences
Sequ_dist_LV <- stringdistmatrix( # Takes 5s
  a = setNames(seqtab.final$Sequence, seqtab.final$ASV_ID),
  b = setNames(Table_16S_279_unique$Sequence, Table_16S_279_unique$ASV.ID),
  method = "lv", useNames = "names")

# Transform to long format
Sequ_dist_LV <- melt(Sequ_dist_LV, varnames = c("ASV_ID", "Ref_species"),
                     value.name = "LV_dist")

################################################################################
#
# Extract relevant Levenshtein taxonomy & calculate evolutionary DNA distance
#
################################################################################

# ------------------------------------------------------------------------------
# Extract the minimum LV distance per ASV to any expected sequence
# ------------------------------------------------------------------------------

# Build dataframe with ASV_ID and first minimum LV distance to any expected seq
Min_sequ_dist <- Sequ_dist_LV %>% 
  group_by(ASV_ID) %>% 
  # Selects the first entry with lowest LV distance, ignoring ties
  slice(which.min(LV_dist)) %>% 
  mutate(ASV_ID = as.character(ASV_ID)) %>%
  rename("ID" = "ASV_ID", "Min_LV_dist" = "LV_dist")

# Extract ambiguous annotations
Min_sequ_dist_ambig <- Sequ_dist_LV %>% 
  group_by(ASV_ID) %>% 
  filter(rank(LV_dist, ties.method = "min") == 1) %>%
  # Remove rows that received the same species annotation from different copies
  distinct() %>%
  mutate(ASV_ID = as.character(ASV_ID),
         Ref_species = as.character(Ref_species)) %>%
  rename("ID" = "ASV_ID", "Min_LV_dist" = "LV_dist") %>%
  as.data.frame()

# Extract sequences with ambiguous annotations and LV <= 8
Ambig_ASVs <- Min_sequ_dist_ambig[duplicated(Min_sequ_dist_ambig$ID) & 
                                    Min_sequ_dist_ambig$Min_LV_dist < 9, "ID"]
# All ambiguous sequences are annotated as E. coli or S. enterica
Min_sequ_dist_ambig[Min_sequ_dist_ambig$ID %in% Ambig_ASVs, ]

# Create ASV table with relative abundances
seqtab.final.rel <- seqtab.final %>%
  mutate(across(S.ID$all, ~ .x/sum(.x)))
# Check maximum relative abundance of ambiguous sequences: 2.3%
rowMax(seqtab.final.rel[seqtab.final.rel$ASV_ID %in% Ambig_ASVs, S.ID$all])
rm(seqtab.final.rel)

# ------------------------------------------------------------------------------
# Evolutionary DNA distance for ambiguous sequences
# ------------------------------------------------------------------------------

# Calculate evolutionary distance
Sequ_dist_evol <- dist.dna(as.DNAbin(strsplit(setNames(c(
  Table_16S_279_unique[grepl("Salm|Esch", 
                             Table_16S_279_unique$ASV.ID), "Sequence"],
  seqtab.final[seqtab.final$ASV_ID %in% Ambig_ASVs, "Sequence"]),
  nm = c(Table_16S_279_unique[grepl("Salm|Esch", 
                                    Table_16S_279_unique$ASV.ID), "ASV.ID"], 
         Ambig_ASVs)), split = "")))
# Transform to matrix and long format
Sequ_dist_evol <- as.matrix(Sequ_dist_evol)
Sequ_dist_evol <- melt(Sequ_dist_evol, varnames = c("ID", "Ref_species"), 
                       value.name = "Evol_dist")
# Filter only distances between ASVs and Ref_species - only clear cases left
Sequ_dist_evol <- Sequ_dist_evol %>% 
  filter(grepl("ASV", ID), !grepl("ASV", Ref_species))
# Extract minimum annotation per ambiguous ASV
Min_sequ_evol_dist <- Sequ_dist_evol %>% 
  group_by(ID) %>% 
  slice(which.min(Evol_dist)) %>% 
  mutate(ID = as.character(ID)) %>%
  rename("Min_evol_dist" = "Evol_dist", "Ref_evol" = "Ref_species")
# Add evol annotation to LV annotation
Min_sequ_dist <- Min_sequ_dist %>% 
  merge(., Min_sequ_evol_dist, by = "ID", all.x = TRUE) %>% 
  mutate(Ref_species = as.character(Ref_species),
         Ref_evol = as.character(Ref_evol))
# Replace LV annotation with evol annotation
Min_sequ_dist$Ref_species[!is.na(Min_sequ_dist$Ref_evol)] <- 
  Min_sequ_dist$Ref_evol[!is.na(Min_sequ_dist$Ref_evol)]
Min_sequ_dist <- Min_sequ_dist[, c("ID", "Ref_species", "Min_LV_dist")]

# ------------------------------------------------------------------------------
# Add exact taxonomy with 0-10 mismatches to seqtab.final
# ------------------------------------------------------------------------------

# Create wide format of the minimum LV distance between obs. and exp. sequences
Min_sequ_dist_wide <- Min_sequ_dist %>% 
  arrange(Min_LV_dist) %>%
  pivot_wider(names_from = Min_LV_dist, values_from = Ref_species)
colnames(Min_sequ_dist_wide) <- c("ASV_ID", paste0(
  "LV_", colnames(Min_sequ_dist_wide)[2:ncol(Min_sequ_dist_wide)]))
seqtab.final <- merge(seqtab.final, Min_sequ_dist_wide[, 1:12], # ASV_ID, LV0-10
  by = "ASV_ID", all = TRUE)

# Fill "empty" annotations with larger LV, that are already given by smaller LV
seqtab.final$LV_1[!is.na(seqtab.final$LV_0)] <- 
  seqtab.final$LV_0[!is.na(seqtab.final$LV_0)]
seqtab.final$LV_2[!is.na(seqtab.final$LV_1)] <- 
  seqtab.final$LV_1[!is.na(seqtab.final$LV_1)]
seqtab.final$LV_3[!is.na(seqtab.final$LV_2)] <- 
  seqtab.final$LV_2[!is.na(seqtab.final$LV_2)]
seqtab.final$LV_4[!is.na(seqtab.final$LV_3)] <- 
  seqtab.final$LV_3[!is.na(seqtab.final$LV_3)]
seqtab.final$LV_5[!is.na(seqtab.final$LV_4)] <- 
  seqtab.final$LV_4[!is.na(seqtab.final$LV_4)]
seqtab.final$LV_6[!is.na(seqtab.final$LV_5)] <- 
  seqtab.final$LV_5[!is.na(seqtab.final$LV_5)]
seqtab.final$LV_7[!is.na(seqtab.final$LV_6)] <- 
  seqtab.final$LV_6[!is.na(seqtab.final$LV_6)]
seqtab.final$LV_8[!is.na(seqtab.final$LV_7)] <- 
  seqtab.final$LV_7[!is.na(seqtab.final$LV_7)]
seqtab.final$LV_9[!is.na(seqtab.final$LV_8)] <- 
  seqtab.final$LV_8[!is.na(seqtab.final$LV_8)]
seqtab.final$LV_10[!is.na(seqtab.final$LV_9)] <- 
  seqtab.final$LV_9[!is.na(seqtab.final$LV_9)]

# ------------------------------------------------------------------------------
# Export seqtab.final of spike-in samples and negative controls with Q buffer
# ------------------------------------------------------------------------------

# Select relevant samples
otus <- seqtab.final[, c("ASV_ID", intersect(
  c(S.ID$spike, S.ID$NEG2), Metafile[Metafile$buffer == "Z", "Sample_ID"]),
  "LV_4")]
# Remove empty ASVs
otus <- otus %>% filter(rowSums(across(where(is.numeric))) > 0)
# Export ASV table
colnames(otus) # for data upload
Metafile[Metafile$Sample_ID %in% colnames(otus), ]
rm(otus)

rm(Ambig_ASVs, Min_sequ_dist_ambig, Min_sequ_dist_wide, Min_sequ_evol_dist, 
   Sequ_dist_evol)
# R objects after running this script:
# "Exp_dist_279", "Exp_dist_279_mean", "Min_sequ_dist", "Sequ_dist_LV", 
#   "Table_16S", "Table_16S_279", "Table_16S_279_unique"
