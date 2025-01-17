################################################################################
#
# Zymo-Qiagen comparison mock analysis
# DADA2 processing of forward FASTQ files
# Based on IEM SOP for microbiome FASTQ processing - denoising with DADA2 (v1.16)
#   based on: https://benjjneb.github.io/dada2/tutorial.html (2021-07-14)
#
################################################################################

# Set project path
path_proj <- 
  "C:/Users/rauerlui/Projects"
setwd(paste0(path_proj, "/Extraction-bias-correction"))

# ------------------------------------------------------------------------------
# Set DADA2 parameters
# ------------------------------------------------------------------------------

# Path of input FASTQ files (must be demultiplexed, can be fastq.gz)
path_input <- paste0(
  path_proj, "/Input/FASTQ/")
# Path for output (Working directory, sequence table, figures etc.)
path_output <- paste0(
  path_proj, "/Output/")
# Folder for filtered FASTQ files
fastq_folder <- "FASTQ_filt_forward_Inf/"
# Set pattern for distinguishing forward and reverse FASTQ file names, e.g. 
#   for SAMPLE-NAME1_R1_001.fastq and SAMPLE-NAME1_R2_001.fastq
forward_pattern <- "R1_001.fastq"

# Determine cutting parameters through QIIME2 visualisation (first value for 
#   forward reads, second value for reverse reads)
truncLen_par <- c(299)
trimLeft_par <- c(20)
# These settings should be fine for our standard ZIEL V1-V3 sequencing.
# The truncLen_par should be changed according to the sequenced region and 
#   quality profile plot (see below), and should leave enough overlap bases 
#   for merging forward and reverse reads (> 20 bp).
# Example for calculating overlap: 
#   E. coli primer positions: 8-27 forward (20 bp), 517-534 reverse (17 bp)
#   Sequencing: 2x300 bp (incl. primers) -> 300+300-(534-8) = 74 bp overlap
#   With current settings: 299+280-(534-8) = 53 bp overlap (= 74-21bp)
# The trimLeft_par cuts the primers. It can be set to zero if primers are 
#   already removed.

# Optional parameters:
truncQ_par <- c(2) # Truncate reads at the first instance of a quality score 
#   less than or equal to truncQ_par. Set lower if you're losing too many reads.
maxEE_par <- c(Inf) # After truncation, reads with higher than maxEE_par
#   "expected errors" will be discarded. Set higher if you're losing to many reads.
maxLen_par <- c(Inf) # Remove reads with length greater than maxLen_par 
#   before trimming and truncation.
minLen_par <- c(20) # Remove reads with length less than minLen_par after 
#   trimming and truncation.
rm.phix_par <- TRUE # Discard reads that match against the phiX genome 
#   (a viral spike-in).

# ------------------------------------------------------------------------------
# Installation & packages
# ------------------------------------------------------------------------------

# Install required packages
if (!"dada2" %in% row.names(installed.packages())) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("dada2", version = "3.11")
}

# Install remaining packages
required_packages <- c("tidyverse", "reshape2")
install.packages(setdiff(required_packages, rownames(installed.packages())))
rm(required_packages)

# Load required packages 
library(dada2)
library(tidyverse)
library(reshape2)
packageVersion("dada2") # 1.16.0

# ------------------------------------------------------------------------------
# Load the FASTQ files
# ------------------------------------------------------------------------------

# List all files in the input directory
list.files(path_input)

# Forward and reverse fastq filenames 
fnFs <- sort(list.files(path_input, pattern = forward_pattern, 
                        full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# ------------------------------------------------------------------------------
# Inspect read quality profiles
# ------------------------------------------------------------------------------

# Interpretation: 
# Grey-scale background: heatmap of each quality score at each base position
# Green line: mean quality score at each position
# Orange line: quartiles of the quality score distribution
# Red line: proportion of reads that extend to at least that position 
#   (not appearing if all sequences are equally long)
# Blue line: current settings for cutting sequences
# Change truncLen_par and trimLeft_par accordingly if needed.

# Quality profiles for forward reads from the first 2 samples
plotQualityProfile(fnFs[c(1:2, 93:94)]) +
  geom_vline(aes(xintercept = trimLeft_par[1]), colour = "blue", alpha = 0.5) +
  geom_vline(aes(xintercept = truncLen_par[1]), colour = "blue", alpha = 0.5)

# ------------------------------------------------------------------------------
# Quality filtering and sequence trimming/truncating
# ------------------------------------------------------------------------------

# Create output directory
if (!dir.exists(paste0(path_output, fastq_folder))){
  dir.create(paste0(path_output, fastq_folder))}

# Create names for filtered files and for a subdirectory of your output folder
filtFs <- file.path(path_output, fastq_folder, 
                    paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

# Quality filtering
out <- filterAndTrim(
  # R object names:
  fwd = fnFs, filt = filtFs,
  # Filtering parameters: 
  truncLen = truncLen_par, trimLeft = trimLeft_par,
  truncQ = truncQ_par, maxEE = maxEE_par,
  maxLen = maxLen_par, minLen = minLen_par, rm.phix = rm.phix_par, 
  # Runtime/output parameters:
  multithread = FALSE, verbose = TRUE)

# ------------------------------------------------------------------------------
# Actual denoising
# ------------------------------------------------------------------------------

# Learn the error rates
errF <- learnErrors(filtFs, multithread = FALSE, verbose = 1, nbases = 1e+09)

# Visualize the estimated error rates
# Interpretation: 
# Black dots: observed error rates by quality score for each base transition.
# Black line: estimated error rates of the DADA2 algorithm
# Red line: expected error rates
# The black line should represent well the black dots.
plotErrors(errF, nominalQ = TRUE)

# Sample inference
dadaFs <- dada(filtFs[file.exists(filtFs)], err = errF, multithread = FALSE)
seqtab <- makeSequenceTable(dadaFs)

# ------------------------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------------------------

# Number of samples & Number of ASVs
dim(seqtab) 
# Distribution of sequence lengths
table(nchar(getSequences(seqtab))) # all sequences have 279 bp

# Track reads through the pipeline (here a little advanced version that works
#   even if no reads in a sample pass the filter criteria)
getN <- function(x) sum(getUniques(x))
track <- list(data.frame(out, row.names = sample.names), 
              data.frame(sapply(dadaFs, getN)))
track <- lapply(track, function(x) rownames_to_column(x, "sample"))
track <- track %>% reduce(merge, by = "sample", all = TRUE)
track <- track %>% column_to_rownames(var = "sample")
colnames(track) <- c("input", "filtered", "denoisedF")
# Mean percentage of recovered reads per step
round(colMeans(track/track$input*100, na.rm = TRUE), 1) # 86.6% of reads survive
round(colMeans(
  track[!grepl("Wasser", rownames(track)), ]/
    track[!grepl("Wasser", rownames(track)), ]$input*100, na.rm = TRUE), 1)
# Excluding water controls, 87.8% of reads survive

# Check singletons (data has samples as rows and ASVs as columns)
table(apply(seqtab, 2, max) == 1) # No singletons present

# Create an ID for each ASV
ASV_ID <- paste0("ASV_", formatC(seq(ncol(seqtab)), 
                                 width = nchar(ncol(seqtab)), flag = "0"))

# ------------------------------------------------------------------------------
# Final data format for next steps
# ------------------------------------------------------------------------------

# Create final ASV table
seqtab.final <- data.frame(ASV_ID = ASV_ID, t(seqtab), check.names = FALSE)
# Move the sequence column to the end
seqtab.final["Sequence"] <- row.names(seqtab.final)
# Remove rownames
row.names(seqtab.final) <- NULL

# ------------------------------------------------------------------------------
# Save R objects
# ------------------------------------------------------------------------------

###save(list = ls(), file = "R_QZmock_DADA2-forw_2024-08-06.RData")
###load("R_QZmock_DADA2-forw_2024-08-06.RData", verbose = TRUE)

rm(dadaFs, errF, fastq_folder, filtFs, fnFs, forward_pattern, getN, maxEE_par, 
   maxLen_par, minLen_par, out, path_input, path_output, rm.phix_par, seqtab, 
   track, trimLeft_par, truncLen_par, truncQ_par)

# R objects after running this script:
# "ASV_ID", "path_proj", "sample.names", "seqtab.final"
