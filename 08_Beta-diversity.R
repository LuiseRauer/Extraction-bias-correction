################################################################################
#
# Zymo-Qiagen comparison mock analysis
# Beta diversity
#
################################################################################

#devtools::install_github("cmartin/ggConvexHull")

# Load required packages
library(tidyverse)
library(vegan) # PCoA
library(metacal) # clr function

################################################################################
#
# Effect of extraction protocols in even mock (plots: BC, p-value: both)
#
################################################################################

# Create an ASV table on species level
seqtab_rel_mock <- seqtab.final %>%
  # Exclude NAs and spikein
  filter(!is.na(LV_4) & !grepl("Allob|Imtech|Truep", LV_4)) %>%
  # Summarise by expected species
  group_by(LV_4) %>%
  summarise_at(.vars = sample.names, .funs = sum) %>%
  # Scale to 100% to account for absence of chimeras and unknowns
  mutate(across(S.ID$all, ~ .x/sum(.x)))

# Run PCoA
PCoA_mock.even <- wcmdscale(vegdist(as.data.frame(t(
  seqtab_rel_mock[, S.ID$mock.even])), method = "bray"), k = 2, eig = FALSE)
colnames(PCoA_mock.even) <- c("PCoA1", "PCoA2")
Metafile_mock.even <- cbind.data.frame(
  Metafile[Metafile$Sample_ID %in% S.ID$mock.even, ], PCoA_mock.even)

# Dilution
svg("Output/Plots/F04A_BC-even-Dil_2023-06-29.svg", 
       width = 4.075, height = 3) 
ggplot(Metafile_mock.even, aes(x = PCoA1, y = PCoA2, col = Dil, fill = Dil)) +
  geom_point(size = 2) +
  scale_colour_manual("Input cells", values = Col_dil) +
  scale_fill_manual("Input cells", values = Col_dil) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Dilution") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.background = element_blank()) +
  annotate("text", x = 0.52, y = -0.09, hjust = 1, vjust = 0, 
           label = "p = 0.75")
dev.off()

# Protocol
svg("Output/Plots/F04B_BC-even-Prot_2023-08-06.svg", 
       width = 3.92, height = 3)
ggplot(Metafile_mock.even, aes(x = PCoA1, y = PCoA2, col = Prot, fill = Prot)) +
  geom_point(size = 2) +
  scale_colour_manual("Protocol", values = Col_prot) +
  scale_fill_manual("Protocol", values = Col_prot) +
  ggConvexHull::geom_convexhull(alpha = 0.2)+ 
  plot_theme + ggtitle("Extraction protocol") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.52, y = -0.09, hjust = 1, vjust = 0, 
           label = "p = 0.002") 
dev.off()

# Kit
svg("Output/Plots/F04C_BC-even-Kit_2023-08-06.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_mock.even, aes(x = PCoA1, y = PCoA2, col = kit, fill = kit)) +
  geom_point(size = 2) +
  scale_colour_manual("Kit", values = Col_prot[9:10]) +
  scale_fill_manual("Kit", values = Col_prot[9:10]) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Extraction kit") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.88, 0.88),
        legend.background = element_rect(fill = NA)) +
  annotate("text", x = 0.52, y = -0.09, hjust = 1, vjust = 0, 
           label = "p = 0.002")
dev.off()

# Lysis condition
svg("Output/Plots/F04D_BC-even-Lysis_2023-06-29.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_mock.even, aes(x = PCoA1, y = PCoA2, col = str_extract(Prot, "T|S"),
                               fill = str_extract(Prot, "T|S"))) +
  geom_point(size = 2) +
  scale_colour_manual("Speed", values = Col_prot[11:12]) +
  scale_fill_manual("Speed", values = Col_prot[11:12]) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Lysis condition") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.88, 0.88),
        legend.background = element_rect(fill = NA)) +
  annotate("text", x = 0.52, y = -0.09, hjust = 1, vjust = 0, 
           label = "p = 0.002")
dev.off()

# Buffer
svg("Output/Plots/F04E_BC-even-Buffer_2023-06-29.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_mock.even, aes(x = PCoA1, y = PCoA2, col = gsub("._._", "", Prot), 
                               fill = gsub("._._", "", Prot))) +
  geom_point(size = 2) +
  scale_colour_manual("Buffer", values = Col_prot[13:14]) +
  scale_fill_manual("Buffer", values = Col_prot[13:14]) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Buffer") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.88, 0.88),
        legend.background = element_rect(fill = NA)) +
  annotate("text", x = 0.52, y = -0.09, hjust = 1, vjust = 0, 
           label = "p = 0.35")
dev.off()

# ------------------------------------------------------------------------------
# PERMANOVA with Bray-Curtis dissimilarity
# ------------------------------------------------------------------------------

# Dilution, p = 0.75
set.seed(1)
adonis2(vegdist(t(seqtab_rel_mock[, S.ID$mock.even]), method = "bray") ~ 
          Metafile_mock.even$Dil, permutations = 500)
# Protocol, p = 0.002
set.seed(1)
adonis2(vegdist(t(seqtab_rel_mock[, S.ID$mock.even]), method = "bray") ~ 
          Metafile_mock.even$Prot, permutations = 500) 
# Kit, p = 0.002
set.seed(1)
adonis2(vegdist(t(seqtab_rel_mock[, S.ID$mock.even]), method = "bray") ~ 
          Metafile_mock.even$kit, permutations = 500)
# Lysis condition, p = 0.002
set.seed(1)
adonis2(vegdist(t(seqtab_rel_mock[, S.ID$mock.even]), method = "bray") ~ 
          Metafile_mock.even$lysis.prot, permutations = 500)
# Buffer, p = 0.35
set.seed(1)
adonis2(vegdist(t(seqtab_rel_mock[, S.ID$mock.even]), method = "bray") ~ 
          Metafile_mock.even$buffer, permutations = 500)

# ------------------------------------------------------------------------------
# PERMANOVA with Aitchison distance
# ------------------------------------------------------------------------------

# Calculate Aitchison distance excluding NAs and spike-in (plot) 
seqtab_CLR_mock <- seqtab.final %>%
  # Exclude NAs and spikein
  filter(!is.na(LV_4) & !grepl("Allob|Imtech|Truep", LV_4)) %>%
  # Summarise by expected species
  group_by(LV_4) %>%
  summarise_at(.vars = sample.names, .funs = sum) %>%
  column_to_rownames("LV_4")

# Replace zeros by 0.5 and apply CLR transformation
seqtab_CLR_mock[seqtab_CLR_mock == 0] <- 0.5
seqtab_CLR_mock <- apply(seqtab_CLR_mock, 2, function(x) clr(x))
round(colSums(seqtab_CLR_mock), 10) # check

# Dilution, p = 0.85
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$mock.even]), method = "euclidean") ~ 
          Metafile_mock.even$Dil, permutations = 500)
# Protocol, p = 0.002
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$mock.even]), method = "euclidean") ~ 
         Metafile_mock.even$Prot, permutations = 500)
# Kit, p = 0.002
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$mock.even]), method = "euclidean") ~ 
          Metafile_mock.even$kit, permutations = 500)
# Lysis condition, p = 0.002
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$mock.even]), method = "euclidean") ~ 
         Metafile_mock.even$lysis.prot, permutations = 500)
# Buffer, p = 0.26
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$mock.even]), method = "euclidean") ~ 
         Metafile_mock.even$buffer, permutations = 500)

################################################################################
#
# Effect of extraction protocols in log mock (plots: BC, p-value: both)
#
################################################################################

# Run PCoA
PCoA_mock.log <- wcmdscale(vegdist(as.data.frame(t(
  seqtab_rel_mock[, S.ID$mock.log])), method = "bray"), k = 2, eig = FALSE)
colnames(PCoA_mock.log) <- c("PCoA1", "PCoA2")
Metafile_mock.log <- cbind.data.frame(
  Metafile[Metafile$Sample_ID %in% S.ID$mock.log, ], PCoA_mock.log)

# Dilution
svg("Output/Plots/S05A_BC-log-Dil_2023-06-29.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_mock.log, aes(x = PCoA1, y = PCoA2, col = Dil, fill = Dil)) +
  geom_point(size = 2) +
  scale_colour_manual("Input cells", values = Col_dil) +
  scale_fill_manual("Input cells", values = Col_dil) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Dilution") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.79, 0.88),
        legend.background = element_rect(fill = NA)) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1)) +
  annotate("text", x = 0.125, y = -0.013, hjust = 1, vjust = 0, 
           label = "p = 0.066") 
dev.off()

# Protocol
svg("Output/Plots/S05B_BC-log-Prot_2023-08-06.svg", 
       width = 3.92, height = 3)
ggplot(Metafile_mock.log, aes(x = PCoA1, y = PCoA2, col = Prot, fill = Prot)) +
  geom_point(size = 2) +
  scale_colour_manual("Protocol", values = Col_prot) +
  scale_fill_manual("Protocol", values = Col_prot) +
  ggConvexHull::geom_convexhull(alpha = 0.2)+ 
  plot_theme + ggtitle("Extraction protocol") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1)) +
  annotate("text", x = 0.125, y = -0.013, hjust = 1, vjust = 0, 
           label = "p = 0.12") 
dev.off()

# Kit
svg("Output/Plots/S05C_BC-log-Kit_2023-08-06.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_mock.log, aes(x = PCoA1, y = PCoA2, col = kit, fill = kit)) +
  geom_point(size = 2) +
  scale_colour_manual("Kit", values = Col_prot[9:10]) +
  scale_fill_manual("Kit", values = Col_prot[9:10]) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Extraction kit") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.88, 0.88),
        legend.background = element_rect(fill = NA)) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1)) +
  annotate("text", x = 0.125, y = -0.013, hjust = 1, vjust = 0, 
           label = "p = 0.34") 
dev.off()

# Lysis condition
svg("Output/Plots/S05D_BC-log-Lysis_2023-06-29.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_mock.log, aes(x = PCoA1, y = PCoA2, col = str_extract(Prot, "T|S"),
                              fill = str_extract(Prot, "T|S"))) +
  geom_point(size = 2) +
  scale_colour_manual("Speed", values = Col_prot[11:12]) +
  scale_fill_manual("Speed", values = Col_prot[11:12]) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Lysis condition") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.88, 0.88),
        legend.background = element_rect(fill = NA)) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1)) +
  annotate("text", x = 0.125, y = -0.013, hjust = 1, vjust = 0, 
           label = "p = 0.008") 
dev.off()

# Buffer
svg("Output/Plots/S05E_BC-log-Buffer_2023-06-29.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_mock.log, aes(x = PCoA1, y = PCoA2, col = gsub("._._", "", Prot), 
                              fill = gsub("._._", "", Prot))) +
  geom_point(size = 2) +
  scale_colour_manual("Buffer", values = Col_prot[13:14]) +
  scale_fill_manual("Buffer", values = Col_prot[13:14]) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Buffer") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.88, 0.88),
        legend.background = element_rect(fill = NA)) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1)) +
  annotate("text", x = 0.125, y = -0.013, hjust = 1, vjust = 0, 
           label = "p = 0.97") 
dev.off()

# ------------------------------------------------------------------------------
# PERMANOVA with Bray-Curtis dissimilarity
# ------------------------------------------------------------------------------

# Dilution, p = 0.066
set.seed(1)
adonis2(vegdist(t(seqtab_rel_mock[, S.ID$mock.log]), method = "bray") ~ 
          Metafile_mock.log$Dil, permutations = 500)
# Protocol, p = 0.12
set.seed(1)
adonis2(vegdist(t(seqtab_rel_mock[, S.ID$mock.log]), method = "bray") ~ 
          Metafile_mock.log$Prot, permutations = 500)
# Kit, p = 0.34
set.seed(1)
adonis2(vegdist(t(seqtab_rel_mock[, S.ID$mock.log]), method = "bray") ~ 
          Metafile_mock.log$kit, permutations = 500)
# Lysis condition, p = 0.008
set.seed(1)
adonis2(vegdist(t(seqtab_rel_mock[, S.ID$mock.log]), method = "bray") ~ 
          Metafile_mock.log$lysis.prot, permutations = 500)
# Buffer, p = 0.97
set.seed(1)
adonis2(vegdist(t(seqtab_rel_mock[, S.ID$mock.log]), method = "bray") ~ 
          Metafile_mock.log$buffer, permutations = 500)

# ------------------------------------------------------------------------------
# PERMANOVA with Aitchison distance
# ------------------------------------------------------------------------------

# Dilution, p = 0.086
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$mock.log]), method = "euclidean") ~ 
          Metafile_mock.log$Dil, permutations = 500)
# Protocol, p = 0.22
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$mock.log]), method = "euclidean") ~ 
         Metafile_mock.log$Prot, permutations = 500)
# Kit, p = 0.50
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$mock.log]), method = "euclidean") ~ 
         Metafile_mock.log$kit, permutations = 500)
# Lysis condition, p = 0.21
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$mock.log]), method = "euclidean") ~ 
         Metafile_mock.log$lysis.prot, permutations = 500)
# Buffer, p = 0.70
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$mock.log]), method = "euclidean") ~ 
         Metafile_mock.log$buffer, permutations = 500)

################################################################################
#
# Beta diversity for spike-in samples
#
################################################################################

# Create an ASV table on species level
seqtab_rel_spike <- seqtab.final %>%
  # Exclude NAs and keep only spikein taxa
  filter(!is.na(LV_4) & grepl("Allob|Imtech|Truep", LV_4)) %>%
  # Summarise by expected species
  group_by(LV_4) %>%
  summarise_at(.vars = sample.names, .funs = sum) %>%
  # Transform to relative abundance
  mutate(across(S.ID$all, ~ .x/sum(.x))) %>%
  # Exclude the one spikein sample that does not contain any spike-in taxa
  select(-`550-BiomPsy-979`)

# Run PCoA
PCoA_spike <- wcmdscale(vegdist(t(
  seqtab_rel_spike[, S.ID$spike[S.ID$spike != "550-BiomPsy-979"]]),
  method = "bray"), k = 2, eig = FALSE)
colnames(PCoA_spike) <- c("PCoA1", "PCoA2")
Metafile_spike <- cbind.data.frame(
  Metafile[Metafile$Sample_ID %in% S.ID$spike[S.ID$spike != "550-BiomPsy-979"], ], 
  PCoA_spike)

# Dilution
svg("Output/Plots/S05F_BC-spike-Dil_2023-06-29.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_spike, aes(x = PCoA1, y = PCoA2, col = Dil, fill = Dil)) +
  geom_point(size = 2) +
  scale_colour_manual("Input cells", values = Col_dil) +
  scale_fill_manual("Input cells", values = Col_dil) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Dilution") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.25, 0.155),
        legend.background = element_rect(fill = NA)) +
  scale_y_continuous(breaks = c(-0.06, -0.03, 0, 0.03)) +
  annotate("text", x = 0.38, y = -0.067, hjust = 1, vjust = 0, 
           label = "p = 0.31") 
dev.off()

# Protocol
svg("Output/Plots/S05G_BC-spike-Prot_2023-08-06.svg", 
       width = 3.92, height = 3)
ggplot(Metafile_spike, aes(x = PCoA1, y = PCoA2, col = Prot, fill = Prot)) +
  geom_point(size = 2) +
  scale_colour_manual("Protocol", values = Col_prot) +
  scale_fill_manual("Protocol", values = Col_prot) +
  ggConvexHull::geom_convexhull(alpha = 0.2)+ 
  plot_theme + ggtitle("Extraction protocol") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = c(-0.06, -0.03, 0, 0.03)) +
  annotate("text", x = 0.38, y = -0.067, hjust = 1, vjust = 0, 
           label = "p = 0.012") 
dev.off()

# Kit
svg("Output/Plots/S05H_BC-spike-Kit_2023-08-06.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_spike, aes(x = PCoA1, y = PCoA2, col = kit, fill = kit)) +
  geom_point(size = 2) +
  scale_colour_manual("Kit", values = Col_prot[9:10]) +
  scale_fill_manual("Kit", values = Col_prot[9:10]) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Extraction kit") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.13, 0.155),
        legend.background = element_rect(fill = NA)) +
  scale_y_continuous(breaks = c(-0.06, -0.03, 0, 0.03)) +
  annotate("text", x = 0.38, y = -0.067, hjust = 1, vjust = 0, 
           label = "p = 0.012") 
dev.off()

# Lysis condition
svg("Output/Plots/S05I_BC-spike-Lysis_2023-06-29.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_spike, aes(x = PCoA1, y = PCoA2, col = str_extract(Prot, "T|S"),
                              fill = str_extract(Prot, "T|S"))) +
  geom_point(size = 2) +
  scale_colour_manual("Speed", values = Col_prot[11:12]) +
  scale_fill_manual("Speed", values = Col_prot[11:12]) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Lysis condition") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.13, 0.155),
        legend.background = element_rect(fill = NA)) +
  scale_y_continuous(breaks = c(-0.06, -0.03, 0, 0.03)) +
  annotate("text", x = 0.38, y = -0.067, hjust = 1, vjust = 0, 
           label = "p = 0.4") 
dev.off()

# Buffer
svg("Output/Plots/S05J_BD-spike-Buffer_2023-06-29.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_spike, aes(x = PCoA1, y = PCoA2, col = gsub("._._", "", Prot), 
                              fill = gsub("._._", "", Prot))) +
  geom_point(size = 2) +
  scale_colour_manual("Buffer", values = Col_prot[13:14]) +
  scale_fill_manual("Buffer", values = Col_prot[13:14]) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Buffer") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.13, 0.155),
        legend.background = element_rect(fill = NA)) +
  scale_y_continuous(breaks = c(-0.06, -0.03, 0, 0.03)) +
  annotate("text", x = 0.38, y = -0.067, hjust = 1, vjust = 0, 
           label = "p = 0.8") 
dev.off()

# ------------------------------------------------------------------------------
# PERMANOVA with Bray-Curtis dissimilarity
# ------------------------------------------------------------------------------

# Dilution, p = 0.31
set.seed(1)
adonis2(vegdist(t(seqtab_rel_spike[, S.ID$spike[S.ID$spike != "550-BiomPsy-979"]]), 
                method = "bray") ~ Metafile_spike$Dil, permutations = 500)
# Protocol, p = 0.012
set.seed(1)
adonis2(vegdist(t(seqtab_rel_spike[, S.ID$spike[S.ID$spike != "550-BiomPsy-979"]]), 
                method = "bray") ~ Metafile_spike$Prot, permutations = 500)
# Kit, p = 0.012
set.seed(1)
adonis2(vegdist(t(seqtab_rel_spike[, S.ID$spike[S.ID$spike != "550-BiomPsy-979"]]), 
                method = "bray") ~ Metafile_spike$kit, permutations = 500)
# Lysis condition, p = 0.4
set.seed(1)
adonis2(vegdist(t(seqtab_rel_spike[, S.ID$spike[S.ID$spike != "550-BiomPsy-979"]]), 
                method = "bray") ~ Metafile_spike$lysis.prot, permutations = 500)
# Buffer, p = 0.8
set.seed(1)
adonis2(vegdist(t(seqtab_rel_spike[, S.ID$spike[S.ID$spike != "550-BiomPsy-979"]]), 
                method = "bray") ~ Metafile_spike$buffer, permutations = 500)

# ------------------------------------------------------------------------------
# PERMANOVA with Aitchison distance
# ------------------------------------------------------------------------------

# Calculate Aitchison distance excluding NAs and spike-in (plot) 
seqtab_CLR_spike <- seqtab.final %>%
  select(LV_4, all_of(S.ID$spike[S.ID$spike != "550-BiomPsy-979"])) %>%
  # Exclude NAs and select spikein taxa
  filter(!is.na(LV_4) & grepl("Allob|Imtech|Truep", LV_4)) %>%
  group_by(LV_4) %>%
  summarise_at(.vars = S.ID$spike[S.ID$spike != "550-BiomPsy-979"], .funs = sum) %>%
  column_to_rownames("LV_4")

# Replace zeros by 0.5
seqtab_CLR_spike[seqtab_CLR_spike == 0] <- 0.5
seqtab_CLR_spike <- apply(seqtab_CLR_spike, 2, function(x) clr(x))
round(colSums(seqtab_CLR_spike), 10) # check

# Dilution, p = 0.68
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_spike[, S.ID$spike[S.ID$spike != "550-BiomPsy-979"]]), 
                method = "euclidean") ~ Metafile_spike$Dil, permutations = 500)
# Protocol, p = 0.028
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_spike[, S.ID$spike[S.ID$spike != "550-BiomPsy-979"]]), 
                method = "euclidean") ~ Metafile_spike$Prot, permutations = 500)
# Kit, p = 0.058
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_spike[, S.ID$spike[S.ID$spike != "550-BiomPsy-979"]]), 
                method = "euclidean") ~ Metafile_spike$kit, permutations = 500)
# Lysis condition, p = 0.13
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_spike[, S.ID$spike[S.ID$spike != "550-BiomPsy-979"]]), 
                method = "euclidean") ~ Metafile_spike$lysis.prot, permutations = 500)
# Buffer, p = 0.43
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_spike[, S.ID$spike[S.ID$spike != "550-BiomPsy-979"]]), 
                method = "euclidean") ~ Metafile_spike$buffer, permutations = 500)

################################################################################
#
# Beta diversity for skin samples with RDP annotation on genus level
#
################################################################################

# Samples 526-BiomPsy-955 & 534-BiomPsy-963 have less reads (protocol Z_T_q)
colSums(seqtab.final[, S.ID$skin])

# Create ASV table with relative abundances
seqtab.final.rel <- seqtab.final %>%
  mutate(across(S.ID$all, ~ .x/sum(.x)))

# Create an ASV table on genus level
seqtab_rel_skin <- seqtab.final.rel %>%
  # Filter 3 Qiagen buffer contaminants & 5 Zymo buffer contaminants
  filter(!grepl("0050|0034|0033|0042|0047|0056|0032|0043", ASV_ID)) %>%
  # Filter genera of buffer contaminants
  filter(!grepl("Methylobacterium|Brucella|Alcaligenes|Paraburkholderia", Genus)) %>%
  filter(!grepl("Nitrospirillum|Herbaspirillum|Aquabacterium", Genus)) %>%
  # Filter spike-in taxa
  filter(!grepl("Truepera|Imtechella|Allobacillus", Genus)) %>%
  # Remove mock taxa from cross-contamination
  filter(is.na(LV_4)) %>%
  # Exclude genera with NA annotation
  filter(Genus != "") %>%
  # Summarise relative abundances on genus level
  group_by(Genus) %>%
  summarise_at(.vars = S.ID$skin, .funs = sum) %>%
  # Scale remaining reads to 100%
  mutate(across(S.ID$skin, ~ .x/sum(.x)))

# Run PCoA
PCoA_skin <- wcmdscale(vegdist(t(seqtab_rel_skin[, S.ID$skin]), 
                               method = "bray"), k = 2, eig = FALSE)
colnames(PCoA_skin) <- c("PCoA1", "PCoA2")
Metafile_skin <- cbind.data.frame(Metafile[Metafile$Sample_ID %in% S.ID$skin, ], 
                                  PCoA_skin)

# Dilution/subject
svg("Output/Plots/S05K_BC-skin-Dil_2023-06-29.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_skin, aes(x = PCoA1, y = PCoA2, col = Dil, fill = Dil)) +
  geom_point(size = 2) +
  scale_colour_manual("Input cells", values = c("grey30", "grey70")) +
  scale_fill_manual("Input cells", values = c("grey30", "grey70")) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Subject") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.22, 0.155),
        legend.background = element_rect(fill = NA)) +
  scale_y_continuous(labels = c("-0.40", "-0.20", "0.00", "0.20")) +
  annotate("text", x = -0.17, y = 0.22, hjust = 0, vjust = 1, 
           label = "p = 0.006") 
dev.off()

# Protocol
svg("Output/Plots/S05L_BC-skin-Prot_2024-09-30.svg", 
       width = 3.92, height = 3)
ggplot(Metafile_skin, aes(x = PCoA1, y = PCoA2, col = Prot, fill = Prot)) +
  geom_point(size = 2) +
  scale_colour_manual("Protocol", values = Col_prot) +
  scale_fill_manual("Protocol", values = Col_prot) +
  ggConvexHull::geom_convexhull(alpha = 0.2)+ 
  plot_theme + ggtitle("Extraction protocol") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = c("-0.40", "-0.20", "0.00", "0.20")) +
  annotate("text", x = -0.17, y = 0.22, hjust = 0, vjust = 1, 
           label = "p = 0.32") 
dev.off()

# Kit
svg("Output/Plots/S05M_BC-skin-Kit_2024-09-30.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_skin, aes(x = PCoA1, y = PCoA2, col = kit, fill = kit)) +
  geom_point(size = 2) +
  scale_colour_manual("Kit", values = Col_prot[9:10]) +
  scale_fill_manual("Kit", values = Col_prot[9:10]) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Extraction kit") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.13, 0.155),
        legend.background = element_rect(fill = NA)) +
  scale_y_continuous(labels = c("-0.40", "-0.20", "0.00", "0.20")) +
  annotate("text", x = -0.17, y = 0.22, hjust = 0, vjust = 1, 
           label = "p = 0.098") 
dev.off()

# Lysis condition
svg("Output/Plots/S05N_BC-skin-Lysis_2024-09-30.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_skin, aes(x = PCoA1, y = PCoA2, col = str_extract(Prot, "T|S"),
                           fill = str_extract(Prot, "T|S"))) +
  geom_point(size = 2) +
  scale_colour_manual("Speed", values = Col_prot[11:12]) +
  scale_fill_manual("Speed", values = Col_prot[11:12]) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Lysis condition") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.13, 0.155),
        legend.background = element_rect(fill = NA)) +
  scale_y_continuous(labels = c("-0.40", "-0.20", "0.00", "0.20")) +
  annotate("text", x = -0.17, y = 0.22, hjust = 0, vjust = 1, 
           label = "p = 0.69") 
dev.off()

# Buffer
svg("Output/Plots/S05O_BC-skin-Buffer_2023-09-30.svg", 
       width = 2.92, height = 3)
ggplot(Metafile_skin, aes(x = PCoA1, y = PCoA2, col = gsub("._._", "", Prot), 
                           fill = gsub("._._", "", Prot))) +
  geom_point(size = 2) +
  scale_colour_manual("Buffer", values = Col_prot[13:14]) +
  scale_fill_manual("Buffer", values = Col_prot[13:14]) +
  ggConvexHull::geom_convexhull(alpha = 0.2) + 
  plot_theme + ggtitle("Buffer") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.13, 0.155),
        legend.background = element_rect(fill = NA)) +
  scale_y_continuous(labels = c("-0.40", "-0.20", "0.00", "0.20")) +
  annotate("text", x = -0.17, y = 0.22, hjust = 0, vjust = 1, 
           label = "p = 0.19") 
dev.off()

# ------------------------------------------------------------------------------
# PERMANOVA with Bray-Curtis dissimilarity
# ------------------------------------------------------------------------------

# Dilution, p = 0.006
set.seed(1)
adonis2(vegdist(t(seqtab_rel_skin[, S.ID$skin]), method = "bray") ~ 
          Metafile_skin$Dil, permutations = 500)
# Protocol, p = 0.32
set.seed(1)
adonis2(vegdist(t(seqtab_rel_skin[, S.ID$skin]), method = "bray") ~ 
          Metafile_skin$Prot, permutations = 500)
# Kit, p = 0.098
set.seed(1)
adonis2(vegdist(t(seqtab_rel_skin[, S.ID$skin]), method = "bray") ~ 
          Metafile_skin$kit, permutations = 500)
# Lysis condition, p = 0.69
set.seed(1)
adonis2(vegdist(t(seqtab_rel_skin[, S.ID$skin]), method = "bray") ~ 
          Metafile_skin$lysis.prot, permutations = 500)
# Buffer, p = 0.19
set.seed(1)
adonis2(vegdist(t(seqtab_rel_skin[, S.ID$skin]), method = "bray") ~ 
          Metafile_skin$buffer, permutations = 500)

# ------------------------------------------------------------------------------
# PERMANOVA with Aitchison distance
# ------------------------------------------------------------------------------

# Calculate Aitchison distance excluding NAs and spike-in (plot) 
seqtab_CLR_skin <- seqtab.final %>%
  # Filter 3 Qiagen buffer contaminants & 5 Zymo buffer contaminants
  filter(!grepl("0050|0034|0033|0042|0047|0056|0032|0043", ASV_ID)) %>%
  # Filter genera of buffer contaminants
  filter(!grepl("Methylobacterium|Brucella|Alcaligenes|Paraburkholderia", Genus)) %>%
  filter(!grepl("Nitrospirillum|Herbaspirillum|Aquabacterium", Genus)) %>%
  # Filter spike-in taxa
  filter(!grepl("Truepera|Imtechella|Allobacillus", Genus)) %>%
  # Remove mock taxa from cross-contamination
  filter(is.na(LV_4)) %>%
  # Exclude genera with NA annotation
  filter(Genus != "") %>%
  # Summarise relative abundances on genus level
  group_by(Genus) %>%
  summarise_at(.vars = sample.names, .funs = sum) %>%
  column_to_rownames("Genus")

# Replace zeros by 0.5
seqtab_CLR_skin[seqtab_CLR_skin == 0] <- 0.5
seqtab_CLR_skin <- apply(seqtab_CLR_skin, 2, function(x) clr(x))
round(colSums(seqtab_CLR_skin), 10) # check

# Dilution, p = 0.18
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$skin]), method = "euclidean") ~ 
          Metafile_skin$Dil, permutations = 500)
# Protocol, p = 1
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$skin]), method = "euclidean") ~ 
         Metafile_skin$Prot, permutations = 500)
# Kit, p = 0.61
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$skin]), method = "euclidean") ~ 
         Metafile_skin$kit, permutations = 500)
# Lysis condition, p = 0.37
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$skin]), method = "euclidean") ~ 
         Metafile_skin$lysis.prot, permutations = 500)
# Buffer, p = 0.92
set.seed(1)
adonis2(vegdist(t(seqtab_CLR_mock[, S.ID$skin]), method = "euclidean") ~ 
         Metafile_skin$buffer, permutations = 500)

rm(PCoA_mock.even, PCoA_mock.log, PCoA_skin, PCoA_spike)

# R objects after running this script:
# "Metafile_mock.even", "Metafile_mock.log", "Metafile_skin", "Metafile_spike",
#   "seqtab_CLR_mock", "seqtab_CLR_skin", "seqtab_CLR_spike", "seqtab_rel_mock", 
#   "seqtab_rel_skin", "seqtab_rel_spike", "seqtab.final.rel
