################################################################################
#
# Zymo-Qiagen comparison mock analysis
# Which extraction protocol is closest to DNA mock?
#
################################################################################

# Load required packages
library(tidyverse)
library(reshape2)
library(vegan)
library(metacal)

################################################################################
#
# Even mock, only considering expected taxa
#
################################################################################

# Define ASV table with expected taxa and count abundances  
seqtab_mock.even <- seqtab.final %>% 
  # Summarise counts on species level
  group_by(LV_4) %>%
  summarise_at(.vars = sample.names, .funs = sum) %>%
  # Select taxonomy information and even cell and DNA mock samples
  select(LV_4, all_of(c(S.ID$mock.even, S.ID$DNA.even))) %>%
  filter(LV_4 %in% Expected_mock$Species) %>%
  column_to_rownames(var = "LV_4")

# ------------------------------------------------------------------------------
# Bray-Curtis dissimilarity
# ------------------------------------------------------------------------------

# Scale read counts to 100%
mock.even_BC <- seqtab_mock.even
mock.even_BC <- apply(mock.even_BC, 2, function(x) x/sum(x)) %>% as.data.frame()
colSums(mock.even_BC) # check sum = 1
# Select samples
mock.even_BC <- mock.even_BC %>% 
  # Calculate compositional mean of the two DNA mock samples
  mutate(DNA = apply(across(S.ID$DNA.even), 1, geom_mean)) %>% 
  # Select only the cell mock samples and the one DNA sample
  select(all_of(c(S.ID$mock.even, "DNA")))
# Calculate BC dissimilarity of all samples, select only distances to DNA sample
mock.even_BC <- data.frame(BC_dist = as.matrix(
  vegdist(t(mock.even_BC), method = "bray"))[, "DNA"]) %>%
  # Remove the distance between DNA and DNA sample
  filter(row.names(.) != "DNA")
# Merge with metafile
mock.even_BC <- merge(mock.even_BC, Metafile, 
                      by.x = 0, by.y = "Sample_ID", all.x = TRUE)

# Plot the result - distance
svg("Output/Plots/S06A_Dist-even_2023-03-26.svg", 
       height = 3.4, width = 4.1)
ggplot(mock.even_BC, aes(x = Prot, y = BC_dist)) +
  geom_boxplot()+ 
  geom_point(size = 2, aes(colour = Dil)) +
  geom_line(aes(group = Dil, colour = Dil), size = 1) +
  scale_colour_manual("Dilution", values = Col_dil) +
  xlab("Protocol") + ylab("Bray-Curtis distance to DNA mock") +
  plot_theme + ggtitle("Even mock") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 0.45)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.4) +
  annotate("text", x = 2.1, y = 0.005, hjust = 1, vjust = 0, label = "No bias",
           size = 3.3)
dev.off()

# Table: distance ranks
mock.even_BC %>% group_by(Prot) %>%
  summarise(BC_dist_median = round(median(BC_dist), 3)) %>%
  mutate(BC_dist_rank = rank(BC_dist_median))

# ------------------------------------------------------------------------------
# Aitchison distance
# ------------------------------------------------------------------------------

# Replace zeros by 0.5
mock.even_Ait <- seqtab_mock.even
mock.even_Ait[mock.even_Ait == 0] <- 0.5 # One Pseudomonas count
# Calculate Aitchison distance
mock.even_Ait <- apply(mock.even_Ait, 2, function(x) clr(x)) %>% as.data.frame()
round(colSums(mock.even_Ait), 10) # check sum = 0
# Select samples
mock.even_Ait <- mock.even_Ait %>% 
  # Calculate mean of 2 DNA mocks 
  # (rowMean on CLR data is identical to compositional mean on count data)
  mutate(DNA = rowMeans(across(S.ID$DNA.even))) %>%
  # Select only the cell mock samples and the one DNA sample
  select(all_of(c(S.ID$mock.even, "DNA")))
# Calculate Aitchison distance of all samples, select only distances to DNA sample
mock.even_Ait <- data.frame(Aitchison_dist = as.matrix(
  vegdist(t(mock.even_Ait), method = "euclidean"))[, "DNA"]) %>%
  # Remove the distance between DNA and DNA sample
  filter(row.names(.) != "DNA")
# Merge with metafile
mock.even_Ait <- merge(mock.even_Ait, Metafile, 
                             by.x = 0, by.y = "Sample_ID", all.x = TRUE)

# Table: distance ranks
mock.even_Ait %>% group_by(Prot) %>%
  summarise(Aitchison_dist_median = round(median(Aitchison_dist), 2)) %>%
  mutate(Aitchison_dist_rank = rank(Aitchison_dist_median))

################################################################################
#
# Staggered mock, only considering expected taxa
#
################################################################################

# Define ASV table with expected taxa and count abundances
seqtab_mock.log <- seqtab.final %>% 
  # Summarise counts on species level
  group_by(LV_4) %>%
  summarise_at(.vars = sample.names, .funs = sum) %>%
  # Select taxonomy information and even cell and DNA mock samples
  select(LV_4, all_of(c(S.ID$mock.log, S.ID$DNA.log))) %>%
  filter(LV_4 %in% Expected_mock$Species) %>%
  # Keep only the 4 species were regularly reads appear (also in DNA):
  filter(grepl("Bac|Pseu|List|Salm", LV_4)) %>%
  column_to_rownames(var = "LV_4")

# ------------------------------------------------------------------------------
# Bray-Curtis
# ------------------------------------------------------------------------------

# Scale read counts to 100%
mock.log_BC <- seqtab_mock.log
mock.log_BC <- apply(mock.log_BC, 2, function(x) x/sum(x)) %>% as.data.frame()
colSums(mock.log_BC) # check sum = 1
# Select samples
mock.log_BC <- mock.log_BC %>% 
  # Calculate compositional mean of the two DNA mock samples
  mutate(DNA = rowMeans(across(S.ID$DNA.log))) %>% 
  # Select only the cell mock samples and the one DNA sample
  select(all_of(c(S.ID$mock.log, "DNA")))
# Calculate BC dissimilarity of all samples, select only distances to DNA sample
mock.log_BC <- data.frame(BC_dist = as.matrix(
  vegdist(t(mock.log_BC), method = "bray"))[, "DNA"]) %>%
  # Remove the distance between DNA and DNA sample
  filter(row.names(.) != "DNA")
# Merge with metafile
mock.log_BC <- merge(mock.log_BC, Metafile, 
                     by.x = 0, by.y = "Sample_ID", all.x = TRUE)

# Plot the result - distance
svg("Output/Plots/S06B_Dist-log_2023-06-28.svg", 
       height = 3.4, width = 2.935)
ggplot(mock.log_BC, aes(x = Prot, y = BC_dist)) +
  geom_boxplot()+ 
  geom_point(size = 2, aes(colour = Dil)) +
  geom_line(aes(group = Dil, colour = Dil), size = 1) +
  scale_colour_manual("Dilution", values = Col_dil) +
  xlab("Protocol") + ylab("Bray-Curtis distance to DNA mock") +
  plot_theme + ggtitle("Staggered mock") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 0.4, 0.1)) +
  coord_cartesian(ylim = c(0, 0.45)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.4) +
  annotate("text", x = 2.1, y = 0.005, hjust = 1, vjust = 0, label = "No bias",
           size = 3.3)
dev.off()

# Table: distance ranks
mock.log_BC %>% group_by(Prot) %>%
  summarise(BC_dist_median = round(median(BC_dist), 3)) %>%
  mutate(BC_dist_rank = rank(BC_dist_median))

# ------------------------------------------------------------------------------
# Aitchison distance
# ------------------------------------------------------------------------------

# Replace zeros by 0.5
mock.log_Ait <- seqtab_mock.log
mock.log_Ait[mock.log_Ait == 0] <- 0.5 # No zeros present
# Calculate Aitchison distance
mock.log_Ait <- apply(mock.log_Ait, 2, function(x) clr(x)) %>% as.data.frame()
round(colSums(mock.log_Ait), 10) # check sum = 0
# Select samples
mock.log_Ait <- mock.log_Ait %>% 
  # Calculate mean of 2 DNA mocks 
  mutate(DNA = rowMeans(across(S.ID$DNA.log))) %>% 
  # Select only the cell mock samples and the one DNA sample
  select(all_of(c(S.ID$mock.log, "DNA")))
# Calculate Aitchison distance of all samples, select only distances to DNA sample
mock.log_Ait <- data.frame(Aitchison_dist = as.matrix(
  vegdist(t(mock.log_Ait), method = "euclidean"))[, "DNA"]) %>%
  # Remove the distance between DNA and DNA sample
  filter(row.names(.) != "DNA")
# Merge with metafile
mock.log_Ait <- merge(mock.log_Ait, Metafile, 
                             by.x = 0, by.y = "Sample_ID", all.x = TRUE)

# Table: distance ranks
mock.log_Ait %>% group_by(Prot) %>%
  summarise(Aitchison_dist_median = round(median(Aitchison_dist), 3)) %>%
  mutate(Aitchison_dist_rank = rank(Aitchison_dist_median))

# R objects after running this script:
# "mock.even_Ait", "mock.even_BC", "mock.log_Ait", "mock.log_BC", 
#   "seqtab_mock.even", "seqtab_mock.log"
