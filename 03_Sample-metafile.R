################################################################################
#
# Zymo-Qiagen comparison mock analysis
# Preparation of sample metadata
#
################################################################################

# Read metafile
Metafile <- read.csv2("Input/Metafile.csv")

# Order variables leves
Metafile$Broad_type <- factor(Metafile$Broad_type, levels = c(
  "Even mock", "Staggered mock", "Spike-in", "Skin", "Neg. control"))
Metafile$Dil <- factor(Metafile$Dil, levels = c(
  "10^8 cells", "10^6 cells", "10^5 cells", "10^4 cells", "D", "6*10^3 cells",
  "Subject 1", "Subject 2", "Pipeline", "A"))

# Define a list of sample IDs for specific tasks
S.ID <- list(
  # all samples
  all = Metafile$Sample_ID,
  # all mock/skin samples excluding NEG1 & NEG2
  sample = Metafile[!grepl("Neg", Metafile$Sample_type), ]$Sample_ID,
  # all NEG2 controls
  NEG2 = Metafile[grepl("Neg.Pipe", Metafile$Sample_type), ]$Sample_ID,
  # all PCR controls
  NEG1 = Metafile[grepl("Neg.PCR", Metafile$Sample_type), ]$Sample_ID,
  # all even cell mock samples
  mock.even = Metafile[Metafile$Sample_type == "Mock.cell.even", ]$Sample_ID,
  # all log cell mock samples 
  mock.log = Metafile[Metafile$Sample_type == "Mock.cell.log", ]$Sample_ID,
  # all mock samples (without skin or control samples)
  mock.all = Metafile[grepl("mock|Spike-in", Metafile$Broad_type), ]$Sample_ID,
  # only skin samples
  skin = Metafile[Metafile$Sample_type %in% c("Skin.P1", "Skin.P2"), ]$Sample_ID,
  # Everything except the skin samples
  notskin = Metafile[!Metafile$Sample_type %in% c("Skin.P1", "Skin.P2"), ]$Sample_ID,
  # only spike samples
  spike = Metafile[Metafile$Sample_type == "Spike", ]$Sample_ID,
  # even/log DNA mock samples
  DNA.even = Metafile[grepl("DNA.even", Metafile$Sample_type), ]$Sample_ID,
  DNA.log = Metafile[grepl("DNA.log", Metafile$Sample_type), ]$Sample_ID
)

# Add ID for the 8 protocols - even mock
for (i in 1:8) { # mock.even.protocol
  S.ID.sub.list <- list(S.ID$mock.even[i+c(0, 8, 16, 24)])
  names(S.ID.sub.list) <- paste0("m.e.p", i)
  S.ID <- append(S.ID, S.ID.sub.list)
}
# Add ID for the 8 protocols - log mock
for (i in 1:8) { # mock.log.protocol
  S.ID.sub.list <- list(S.ID$mock.log[i+c(0, 8)])
  names(S.ID.sub.list) <- paste0("m.l.p", i)
  S.ID <- append(S.ID, S.ID.sub.list)
}
# Add ID for the 8 protocols - spike
for (i in 1:8) { # mock.log.protocol
  S.ID.sub.list <- list(S.ID$spike[i+c(0, 8)])
  names(S.ID.sub.list) <- paste0("m.s.p", i)
  S.ID <- append(S.ID, S.ID.sub.list)
}
# Add ID for the 4 dilutions
for (i in 1:4) { # mock.even.dilution
  S.ID.sub.list <- list(S.ID$mock.even[c(1:8)+(i-1)*8])
  names(S.ID.sub.list) <- paste0("m.e.dil", i)
  S.ID <- append(S.ID, S.ID.sub.list)
}
# Add ID for the 2 dilutions
for (i in 1:2) { # mock.log.dilution
  S.ID.sub.list <- list(S.ID$mock.log[c(1:8)+(i-1)*8])
  names(S.ID.sub.list) <- paste0("m.l.dil", i)
  S.ID <- append(S.ID, S.ID.sub.list)
}
rm(i, S.ID.sub.list)

# R objects after running this script:
# "S.ID", "Metafile"
