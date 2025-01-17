################################################################################
#
# Zymo-Qiagen comparison mock analysis
# Analyse bias per protocol and species using metacal
#
################################################################################

library(tidyverse)
library(metacal) # estimate_bias
library(rstatix) # significance bars
library(ggpubr) # stat_pvalue_manual, as_ggplot
library(cowplot) # plot_grid, get_legend 
library(reshape2) # melt
library(vegan) # vegdist
#install.packages("ggsignif") 
#library(ggsignif)

###save(list = ls(), file = "R_QZmock_scripts1-9_2024-08-28.RData")
###load("R_QZmock_scripts1-9_2024-08-28.RData", verbose = TRUE)

################################################################################
#
# Estimate extraction bias per sample
#
################################################################################

# ------------------------------------------------------------------------------
# Even mock
# ------------------------------------------------------------------------------

# Define "observed" data (cell mocks) & remove sample with 0 Pseudomonas counts
Obs_mock.even <- seqtab_mock.even %>% 
  select(-all_of(S.ID$DNA.even))

# Define "expected" data (mean composition of DNA mocks)
Exp_DNA.even <- data.frame(matrix(
  # Calculate the geometric mean of the two internal DNA mocks
  apply(seqtab_mock.even %>% 
          select(all_of(S.ID$DNA.even)) %>%
          arrange(match(row.names(.), Expected_mock$Species)), 1, geom_mean), 
  nrow = 8, ncol = length(S.ID$mock.even), 
  dimnames = list(Expected_mock$Species, S.ID$mock.even)), 
  check.names = FALSE) %>%
  arrange(row.names(.))

# Check that Observed and Expected data are in same format
identical(row.names(Obs_mock.even), row.names(Exp_DNA.even))
identical(colnames(Obs_mock.even), colnames(Exp_DNA.even))

# Estimate bias in cell mocks per sample
Mock.even_prot_bias <- list()
for (i in c(1:32)) { # 8 protocols * 4 dilutions
  S.ID.metacal <- S.ID$mock.even[i]
  if(i == 31) { # Exclude Pseudomonas
    set.seed(1)
    Mock.even_prot_bias[[i]] <- estimate_bias(
      as.matrix(Obs_mock.even[c(1:5, 7:8), S.ID.metacal, drop = F]), 
      as.matrix(Exp_DNA.even[c(1:5, 7:8), S.ID.metacal, drop = F]), 2, boot = TRUE)
  } else {
    set.seed(1)
    Mock.even_prot_bias[[i]] <- estimate_bias(
      as.matrix(Obs_mock.even[, S.ID.metacal, drop = F]), 
      as.matrix(Exp_DNA.even[, S.ID.metacal, drop = F]), 2, boot = TRUE)
  }
}
# Extract coefficients
Mock.even_prot_bias <- lapply(Mock.even_prot_bias[c(1:32)], function(x) 
  summary(x)$coefficients)
names(Mock.even_prot_bias) <- 
  Metafile[Metafile$Sample_ID %in% S.ID$mock.even, "Sample_ID"]
Mock.even_prot_bias <- Mock.even_prot_bias %>% bind_rows(.id = "Sample_ID")

# ------------------------------------------------------------------------------
# Staggered mock
# ------------------------------------------------------------------------------

# Define "observed" data (cell mocks), no zero counts present
Obs_mock.log <- seqtab_mock.log %>% 
  select(-all_of(S.ID$DNA.log))

# Define "expected" data (mean composition of DNA mocks)
Exp_DNA.log <- data.frame(matrix(
  # Calculate the geometric mean of the two DNA mocks
  apply(seqtab_mock.log %>% 
          select(all_of(S.ID$DNA.log)), 1, geom_mean), 
  nrow = 4, ncol = length(S.ID$mock.log), 
  dimnames = list(row.names(seqtab_mock.log), S.ID$mock.log)), 
  check.names = FALSE) %>%
  arrange(row.names(.))

# Estimate bias in cell mocks
Mock.log_prot_bias <- list()
for (i in 1:16) {
  S.ID.metacal <- S.ID$mock.log[i]
  set.seed(1)
  Mock.log_prot_bias[[i]] <- estimate_bias(
    as.matrix(Obs_mock.log[, S.ID.metacal, drop = F]), 
    as.matrix(Exp_DNA.log[, S.ID.metacal, drop = F]), 2, boot = TRUE)
}

# Extract coefficients
Mock.log_prot_bias <- lapply(Mock.log_prot_bias, function(x) 
  summary(x)$coefficients)
names(Mock.log_prot_bias) <- 
  Metafile[Metafile$Sample_ID %in% S.ID$mock.log, "Sample_ID"]
Mock.log_prot_bias <- Mock.log_prot_bias %>% bind_rows(.id = "Sample_ID")

# ------------------------------------------------------------------------------
# Spike-in mock
# ------------------------------------------------------------------------------

# Define observed data (cell and DNA)
Obs_spike <- seqtab.final %>%
  # Summarise taxonomy at LV_4
  group_by(LV_4) %>%
  summarise_at(.vars = sample.names, .funs = sum) %>%
  # Remove NA annotation
  filter(LV_4 %in% Expected_spike$Species) %>%
  column_to_rownames(var = "LV_4") %>%
  # Select all even mock samples (DNA & cell)
  # EXCLUDE THE SAME SAMPLES AS IN LOWER LINE
  select(all_of(S.ID$spike[S.ID$spike != "550-BiomPsy-979"]))
# Replace zeros by 0.5 (0 for Pseudomonas in 502-BiomPsy-931)
#Obs_spike[Obs_spike == 0] <- 0.5

# Define expected data
Exp_spike <- data.frame(matrix(
  Expected_spike$Abund_16S, 
  nrow = 3, ncol = length(S.ID$spike[!S.ID$spike %in% c("550-BiomPsy-979", "542-BiomPsy-971")]), 
  dimnames = list(Expected_spike$Species, S.ID$spike[!S.ID$spike %in% c("550-BiomPsy-979", "542-BiomPsy-971")])), 
  check.names = FALSE) %>%
  arrange(row.names(.))

# Estimate bias in cell mocks
Spike_prot_bias <- list()
for (i in c(1:14)) {
  S.ID.metacal <- S.ID$spike[c(1:6, 8:14, 16)][i]
  if(S.ID.metacal == "549-BiomPsy-978") {
    set.seed(1)
    Spike_prot_bias[[i]] <- estimate_bias(
      as.matrix(Obs_spike[2:3, S.ID.metacal, drop = F]), 
      as.matrix(Exp_spike[2:3, S.ID.metacal, drop = F]), 2, boot = TRUE)
  } else {
    set.seed(1)
    Spike_prot_bias[[i]] <- estimate_bias(
      as.matrix(Obs_spike[, S.ID.metacal, drop = F]), 
      as.matrix(Exp_spike[, S.ID.metacal, drop = F]), 2, boot = TRUE)
  }
}

# Extract coefficients
Spike_prot_bias <- lapply(Spike_prot_bias, function(x) 
  summary(x)$coefficients)
names(Spike_prot_bias) <- 
  Metafile[Metafile$Sample_ID %in% S.ID$spike[c(1:6, 8:14, 16)], "Sample_ID"]
Spike_prot_bias <- Spike_prot_bias %>% bind_rows(.id = "Sample_ID")

# ------------------------------------------------------------------------------
# Combine the data
# ------------------------------------------------------------------------------

# Combine the data & define morphology groups
Prot_bias <- rbind.data.frame(
  Mock.even_prot_bias %>% mutate(Mock = "Even mock"), 
  Mock.log_prot_bias %>% mutate(Mock = "Staggered mock"), 
  Spike_prot_bias %>% mutate(Mock = "Spike-in mock")) %>% 
  mutate(Morph = case_when(
    grepl("Pseudo|Salmon|Esche|Imtech", taxon) ~ "rods, gram neg.",
    grepl("Bacillus|Lactoba|Allob", taxon, ignore.case = F) ~ "rods, gram pos.",
    grepl("List|Enteroc|Staph|Truep", taxon) ~ "coccoid, gram pos.")) 
# Merge with metafile
Prot_bias <- Prot_bias %>% 
  merge(., Metafile[, c("Sample_ID", "Prot", "Dil")], by = "Sample_ID", 
        all.x = TRUE)
colnames(Prot_bias)[colnames(Prot_bias) == "Prot"] <- "Protocol"
# Extract protocol information
Prot_bias <- Prot_bias %>% mutate(
  Kit = substring(Protocol, 1, 1),
  Lysis = substring(Protocol, 3, 3),
  Buffer = substring(Protocol, 5, 5))
# Extract Gram & shape, combine protocols with different buffers
Prot_bias <- Prot_bias %>% mutate(
  Gram = ifelse(grepl("pos", Morph), "pos", "neg"),
  Shape = ifelse(grepl("rod", Morph), "rod", "coccus"),
  Group = case_when(grepl("Q_S", Protocol) ~ "Q_S_q, Q_S_z",
                    grepl("Q_T", Protocol) ~ "Q_T_q, Q_T_z",
                    grepl("Z_S", Protocol) ~ "Z_S_q, Z_S_z",
                    TRUE ~ "Z_T_q, Z_T_z"))
# Make dilution information numeric 
Prot_bias$Dil <- as.character(Prot_bias$Dil)
Prot_bias$Dil <- gsub(" cells", "", Prot_bias$Dil)
Prot_bias$Dil <- gsub("0\\^", "e", Prot_bias$Dil)
Prot_bias$Dil <- gsub("6\\*1", "6", Prot_bias$Dil)
Prot_bias$Dil <- as.numeric(Prot_bias$Dil)

# Add significance bars
Plot_signif_bars <- Prot_bias %>%
  group_by(Group) %>%
  wilcox_test(estimate ~ Morph, data = .) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>% 
  add_xy_position(x = c("Morph"), fun = "max")
# Fix position of bars in the plot
Plot_signif_bars <- Plot_signif_bars %>%
  mutate(xmin = case_when(xmin == 1 ~ 3,
                          xmin == 2 ~ 8.5),
         xmax = case_when(xmin == 3 & xmax == 2 ~ 14.5,
                          xmin == 3 & xmax == 3 ~ 8.5,
                          TRUE ~ 14.5),
         y.position = case_when(xmin == 3 & xmax == 8.5 ~ 0.9,
                                xmin == 3 & xmax == 14.5 ~ 1.07,
                                xmin == 8.5 & xmax == 14.5 ~ 1.24))

# Plot bias per species with significance bars
svg("Output/Plots/F05_Bias_all_2024-09-30.svg", width = 6, height = 7.5)
set.seed(1)
Prot_bias %>% 
  # Add fake space between morphology groups on x axis
  add_row(taxon = c(" ", "  "), Group = "Q_T_q, Q_T_z", Protocol = "Q_T_q",
          Mock = "Even mock") %>%
  # Fix order of morphology groups, mocks, and taxa
  mutate(Morph = factor(Morph, levels = c("coccoid, gram pos.", "rods, gram pos.",
                                          "rods, gram neg.")),
         Mock = factor(Mock, levels = c("Even mock", "Staggered mock",
                                        "Spike-in mock")),
         taxon2 = paste0(taxon, ".", Mock),
         taxon3 = case_when(
           taxon2 == "Enterococcus_faecalis.Even mock" ~ 1,
           taxon2 == "Listeria_monocytogenes.Even mock" ~ 2,
           taxon2 == "Listeria_monocytogenes.Staggered mock" ~ 3,
           taxon2 == "Staphylococcus_aureus.Even mock" ~ 4,
           taxon2 == "Truepera_radiovictrix.Spike-in mock" ~ 5,
           taxon2 == " .Even mock" ~ 6,
           taxon2 == "Allobacillus_halotolerans.Spike-in mock" ~ 7,
           taxon2 == "Bacillus_subtilis.Even mock" ~ 8,
           taxon2 == "Bacillus_subtilis.Staggered mock" ~ 9,
           taxon2 == "Lactobacillus_fermentum.Even mock" ~ 10,
           taxon2 == "  .Even mock" ~ 11,
           taxon2 == "Escherichia_coli.Even mock" ~ 12,
           taxon2 == "Imtechella_halotolerans.Spike-in mock" ~ 13,
           taxon2 == "Pseudomonas_aeruginosa.Even mock" ~ 14,
           taxon2 == "Pseudomonas_aeruginosa.Staggered mock" ~ 15,
           taxon2 == "Salmonella_enterica.Even mock" ~ 16,
           taxon2 == "Salmonella_enterica.Staggered mock" ~ 17)) %>%
  ggplot(.) +
  geom_hline(yintercept = 1, colour = "grey50") +
  geom_jitter(aes(x = taxon3, y = estimate, fill = taxon, shape = Mock, 
                  size = log10(Dil)), height = 0.2, width = 0.2, col = "black") + 
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10), 
                labels = c("0.01", "0.1", "1", "10")) +
  facet_grid(Group ~ .) +
  coord_cartesian(ylim = c(0.1, 18.5)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = c(3, 8.5, 14.5), 
                     labels = c("Ovoid, gram pos.", "Rods, gram pos.", 
                                "Rods, gram neg.")) +
  scale_fill_manual("Protocol\n\n", values = Col_spec) +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.margin = margin(0.5, 0, 0.5, 0, unit = "cm"),
        axis.text.x = element_text(colour = "black", size = rel(1))) +
  plot_theme + ylab("Extraction bias") +
  scale_shape_manual("Mock community", values = c(21, 22, 24)) +
  scale_size("\n\n\nInput cells [log10]", range = c(1, 4)) +
  guides(fill = guide_legend(order = 1, title.position = "top", ncol = 1,
                             byrow = TRUE,
                             override.aes = list(size = 3, shape = 21)),
  size = guide_legend(order = 2, title.position = "top", nrow = 1,
                      label.position = "bottom", reverse = TRUE),
  shape = guide_legend(order = 3, title.position = "top", nrow = 3,
                       override.aes = list(size = 3)))+
  stat_pvalue_manual(Plot_signif_bars)
dev.off()

# ------------------------------------------------------------------------------
# Legend
# ------------------------------------------------------------------------------

for (i in c("rods, gram neg.", "rods, gram pos.", "coccoid, gram pos.")) {
  p <- Prot_bias %>% 
    mutate(taxon = factor(taxon, levels = c(
      "Enterococcus_faecalis", "Listeria_monocytogenes", "Staphylococcus_aureus", 
      "Truepera_radiovictrix", " ", "Allobacillus_halotolerans", 
      "Bacillus_subtilis", "Lactobacillus_fermentum", "  ", "Escherichia_coli",
      "Imtechella_halotolerans", "Pseudomonas_aeruginosa", "Salmonella_enterica"))) %>%
    filter(Morph == i) %>%
    ggplot() +
    geom_jitter(aes(x = taxon, y = estimate, fill = taxon), colour = "black") + 
    plot_theme +
    theme(legend.text.align = 0,
          legend.margin = margin(0.15, 0, 0.5, 0, unit = "cm"),
          legend.justification = c(0, 1)) +
    guides(fill = guide_legend(order = 1, ncol = 1, byrow = TRUE, 
                               override.aes = list(size = 3, shape = 21)))
  if(i == "rods, gram neg.") {p1 <- p + 
    scale_fill_manual("Rods, gram negative", values = Col_spec, 
                      labels = species_labels)}
  if(i == "rods, gram pos.") {p2 <- p +
    scale_fill_manual("Rods, gram positive", values = Col_spec, 
                      labels = species_labels)}
  if(i == "coccoid, gram pos.") {p3 <- p +
    scale_fill_manual("Ovoid, gram positive", values = Col_spec, 
                      labels = species_labels)}
}
svg("Output/Plots/F05Legend_2024-09-30.svg", width = 1.9, height = 4.5)
plot_grid(as_ggplot(get_legend(p3)), as_ggplot(get_legend(p2)),
          as_ggplot(get_legend(p1)), 
          ncol = 1, align = "hv", rel_heights = c(1.18, 1, 1.15))
dev.off()

rm(p, p1, p2, p3, i)

# ------------------------------------------------------------------------------
# Mann Whitney U test
# ------------------------------------------------------------------------------

# Test for any difference between the three groups
kruskal.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Q_S_q, Q_S_z"))
kruskal.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Z_S_q, Z_S_z"))
kruskal.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Q_T_q, Q_T_z"))
kruskal.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Z_T_q, Z_T_z"))

# Pairwise test for differences
wilcox.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Q_S_q, Q_S_z",
                                            grepl("cocc|pos", Morph)))
wilcox.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Q_S_q, Q_S_z",
                                            grepl("cocc|neg", Morph)))
wilcox.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Q_S_q, Q_S_z",
                                            !grepl("cocc", Morph)))

wilcox.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Z_S_q, Z_S_z",
                                            grepl("cocc|pos", Morph)))
wilcox.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Z_S_q, Z_S_z",
                                            grepl("cocc|neg", Morph)))
wilcox.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Z_S_q, Z_S_z",
                                            !grepl("cocc", Morph)))

wilcox.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Q_T_q, Q_T_z",
                                            grepl("cocc|pos", Morph)))
wilcox.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Q_T_q, Q_T_z",
                                            grepl("cocc|neg", Morph)))
wilcox.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Q_T_q, Q_T_z",
                                            !grepl("cocc", Morph)))

wilcox.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Z_T_q, Z_T_z",
                                            grepl("cocc|pos", Morph)))
wilcox.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Z_T_q, Z_T_z",
                                            grepl("cocc|neg", Morph)))
wilcox.test(estimate ~ Morph, data = filter(Prot_bias, Group == "Z_T_q, Z_T_z",
                                            !grepl("cocc", Morph)))

################################################################################
#
# Correct extraction bias with bias from 10^6 even cell sample
#
################################################################################

# ------------------------------------------------------------------------------
# Generate 8*3 calibration values from previous part of same script
# ------------------------------------------------------------------------------

# Extract bias of the 8 10^6 even cell mock samples
Bias_even_106 <- Mock.even_prot_bias %>% 
  filter(Sample_ID %in% S.ID$m.e.dil2) %>%
  merge(., Metafile[, c("Sample_ID", "Prot")]) %>%
  # Define geometric mean of bias per morphology group
  mutate(Group = case_when(
    grepl("Pseudo|Salmon|Esche", taxon) ~ "rods, gram neg.",
    grepl("Bacillus|Lactoba", taxon, ignore.case = F) ~ "rods, gram pos.",
    grepl("List|Enteroc|Staph", taxon) ~ "coccoid, gram pos.")) %>% 
  group_by(Group, Prot) %>%
  summarise(estimate_gm = geom_mean(estimate))

# ------------------------------------------------------------------------------
# Calibrate other samples
# ------------------------------------------------------------------------------

# Add morphology groups to observed data
Mock.even_cal <- seqtab_mock.even %>% # (cell & DNA)
  mutate(DNA = apply(across(S.ID$DNA.even), 1, geom_mean)) %>% 
  select(all_of(S.ID$mock.even), "DNA") %>%
  rownames_to_column(var = "taxon") %>% melt(., variable.name = "Sample_ID") %>% 
  mutate(Group = case_when(
    grepl("Pseudo|Salmon|Esche", taxon) ~ "rods, gram neg.",
    grepl("Bacillus|Lactoba", taxon, ignore.case = F) ~ "rods, gram pos.",
    grepl("List|Enteroc|Staph", taxon) ~ "coccoid, gram pos."))

# Merge with metafile (add protocol information)
Mock.even_cal <- merge(Mock.even_cal, Metafile[, c("Sample_ID", "Prot")], 
                            by = "Sample_ID", all.x = TRUE)

# Merge bias (3 groups * 8 prot) with obs data ((8 prots * 4 dilutions+1) * 8 species)
Mock.even_cal <- Mock.even_cal %>%
  merge(., Bias_even_106, by = c("Prot", "Group"), all.x = TRUE) 

# Calibrate counts, make relative abundances
Mock.even_cal <- Mock.even_cal %>% 
  # Calibrate observed counts with bias (divide observed values by bias)
  mutate(val_cal = value/estimate_gm) %>%
  # Scale each sample to 100%
  group_by(Sample_ID) %>%
  mutate(val_scale = value/sum(value),
         val_cal_scale = val_cal/sum(val_cal)) %>% #View
  mutate(BD_BC = case_when(is.na(estimate_gm) ~ val_scale, TRUE ~ val_cal_scale), # scaled data for BC
         BD_Ait = case_when(is.na(estimate_gm) ~ value, TRUE ~ val_cal)) # count data for AIT

# ------------------------------------------------------------------------------
# BC
# ------------------------------------------------------------------------------

# Data in final format: Calibrated mock samples (calibrierte rel. ab. in wide format)
Mock.even_cal_BC <- Mock.even_cal %>% 
  select(taxon, Sample_ID, BD_BC) %>%
  pivot_wider(names_from = taxon, values_from = BD_BC) %>%
  column_to_rownames("Sample_ID") %>% 
  t(.) %>% as.data.frame()

# Calculate BC distance, select only column "Expected"
Mock.even_cal_BC_dist <- data.frame(BC_cal = as.matrix(
  vegdist(t(Mock.even_cal_BC), method = "bray"))[, "DNA"]) 
# Correct the distance with the zero Pseudomonas count
Mock.even_cal_BC_dist[row.names(Mock.even_cal_BC_dist) == "502-BiomPsy-931", ] <-
  vegdist(t(Mock.even_cal_BC[!grepl(
    "Pseudo", row.names(Mock.even_cal_BC)), c("502-BiomPsy-931", "DNA")]))
# Remove the distance between DNA and DNA sample
Mock.even_cal_BC_dist <- Mock.even_cal_BC_dist %>%
  filter(row.names(.) != "DNA")
Mock.even_cal_BC_dist <- merge(
  Mock.even_cal_BC_dist, mock.even_BC, by.x = 0, by.y = "Row.names", all = TRUE) %>%
  rename("Original" = BC_dist, "Corrected" = BC_cal)

# Plot the result - boxplot, 10^6, main figure
annotations <- data.frame(
  x = c(0.5), y = c(0.0605),
  Dil = c("10^6 cells"),
  label = c("p = 0.031"))
svg("Output/Plots/F06A_Dist-even-106_2023-08-16.svg", 
       height = 2.24, width = 3) 
Mock.even_cal_BC_dist %>%
  filter(Dil == "10^6 cells") %>%
  select(Dil, Original, Corrected, Prot) %>% melt(.) %>% 
  filter(!grepl("Z_T", Prot)) %>%
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot() + geom_point() +
  geom_line(aes(group = Prot, colour = Prot), size = 2) +
  scale_y_log10() + ylab("BC distance to DNA mock") +
  facet_grid(Dil ~., scales = "free_y") +
  scale_colour_manual(values = Col_prot) + plot_theme + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = rel(1.2), colour = "black")) +
  geom_text(data = annotations, aes(x = x, y = y, label = label), 
            hjust = 0, vjust = 0) 
dev.off()

# Plot the result - boxplot, 10^6, supp figure
annotations <- data.frame(
  x = c(0.5), y = c(0.0605),
  Dil = c("10^6 cells"),
  label = c("p = 0.0078"))
svg("Output/Plots/S07A_Dist-even-106_2023-08-16.svg", 
       height = 2.24, width = 3)
Mock.even_cal_BC_dist %>%
  filter(Dil == "10^6 cells") %>%
  select(Dil, Original, Corrected, Prot) %>% melt(.) %>% 
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot() + geom_point() +
  geom_line(aes(group = Prot, colour = Prot), size = 2) +
  scale_y_log10() + ylab("BC distance to DNA mock") +
  facet_grid(Dil ~., scales = "free_y") +
  scale_colour_manual(values = Col_prot) + plot_theme + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = rel(1.2), colour = "black")) +
  geom_text(data = annotations, aes(x = x, y = y, label = label), 
            hjust = 0, vjust = 0) 
dev.off()

# Plot the result - boxplot, 10^8/5/4, main figure
annotations <- data.frame(
  x = c(0.5, 0.5, 0.5),
  y = c(0.074, 0.0518, 0.0748),
  Dil = factor(c("10^8 cells", "10^5 cells", "10^4 cells"), 
               levels = c("10^8 cells", "10^5 cells", "10^4 cells")),
  label = c("p = 0.031", "p = 0.031", "p = 0.063"))
svg("Output/Plots/F06B_Dist-even-104-8_2023-08-16.svg", 
       height = 6.16, width = 3) # 10^4, 10^5, 10^8
Mock.even_cal_BC_dist %>%
  filter(Dil != "10^6 cells") %>%
  select(Dil, Original, Corrected, Prot) %>% melt(.) %>% 
  filter(!grepl("Z_T", Prot)) %>%
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot() + geom_point() +
  geom_line(aes(group = Prot, colour = Prot), size = 2) +
  scale_y_log10() + ylab("BC distance to DNA mock") +
  facet_grid(Dil ~., scales = "free_y") +
  scale_colour_manual(values = Col_prot) + plot_theme + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = rel(1.2), colour = "black")) +
  geom_text(data = annotations, aes(x = x, y = y, label = label), 
            hjust = 0, vjust = 0)
dev.off()

# Plot the result - boxplot, 10^8/5/4, supp figure
annotations <- data.frame(
  x = c(0.5, 0.5, 0.5),
  y = c(0.074, 0.0518, 0.0748),
  Dil = factor(c("10^8 cells", "10^5 cells", "10^4 cells"), 
               levels = c("10^8 cells", "10^5 cells", "10^4 cells")),
  label = c("p = 0.016", "p = 0.20", "p = 0.016"))
svg("Output/Plots/S07B_Dist-even-104-8_2023-08-16.svg", 
       height = 6.16, width = 3) # 10^4, 10^5, 10^8
Mock.even_cal_BC_dist %>%
  filter(Dil != "10^6 cells") %>%
  select(Dil, Original, Corrected, Prot) %>% melt(.) %>% 
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot() + geom_point() +
  geom_line(aes(group = Prot, colour = Prot), size = 2) +
  scale_y_log10() + ylab("BC distance to DNA mock") +
  facet_grid(Dil ~., scales = "free_y") +
  scale_colour_manual(values = Col_prot) + plot_theme + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = rel(1.2), colour = "black")) +
  geom_text(data = annotations, aes(x = x, y = y, label = label), 
            hjust = 0, vjust = 0)
dev.off()

# ------------------------------------------------------------------------------
# Statistical test
# ------------------------------------------------------------------------------
Even_BC_test <- Mock.even_cal_BC_dist %>%
  select(Dil, Original, Corrected, Prot) %>% melt(.)

# Main, 10^6, p = 0.031
wilcox.test(value ~ variable, data = filter(Even_BC_test, Dil == "10^6 cells",
                                            !grepl("Z_T", Prot)), paired = TRUE)
# Main, 10^8, p = 0.031
wilcox.test(value ~ variable, data = filter(Even_BC_test, Dil == "10^8 cells",
                                            !grepl("Z_T", Prot)), paired = TRUE)
# Main, 10^5, p = 0.031
wilcox.test(value ~ variable, data = filter(Even_BC_test, Dil == "10^5 cells",
                                            !grepl("Z_T", Prot)), paired = TRUE)
# Main, 10^4, p = 0.063
wilcox.test(value ~ variable, data = filter(Even_BC_test, Dil == "10^4 cells",
                                            !grepl("Z_T", Prot)), paired = TRUE)

# Supp, 10^6, p = 0.0078
wilcox.test(value ~ variable, data = filter(Even_BC_test, Dil == "10^6 cells"), 
            paired = TRUE)
# Supp, 10^8, p = 0.016
wilcox.test(value ~ variable, data = Even_BC_test %>% filter(Dil == "10^8 cells"),
            paired = TRUE)
# Supp, 10^5, p = 0.20
wilcox.test(value ~ variable, data = filter(Even_BC_test, Dil == "10^5 cells"), 
            paired = TRUE)
# Supp, 10^4, p = 0.016
wilcox.test(value ~ variable, data = filter(Even_BC_test, Dil == "10^4 cells"), 
            paired = TRUE)

# ------------------------------------------------------------------------------
# Aitchison
# ------------------------------------------------------------------------------

# Calculate Aitchison distance
Mock.even_cal_Ait <- Mock.even_cal %>%
  select(Sample_ID, taxon, BD_Ait) %>%
  pivot_wider(names_from = taxon, values_from = BD_Ait) %>%
  column_to_rownames("Sample_ID") %>% 
  t(.) %>% as.data.frame()
Mock.even_cal_Ait[Mock.even_cal_Ait == 0] <- 0.5
# Aitchison distance
Mock.even_cal_Ait <- apply(Mock.even_cal_Ait, 2, function(x) clr(x))
round(colSums(Mock.even_cal_Ait), 10) # check
# Calculate Aitchison distance, select only column "Expected"
Mock.even_cal_Ait_dist <- data.frame(Aitchison_dist = as.matrix(
  vegdist(t(Mock.even_cal_Ait), method = "euclidean"))[, "DNA"])
# Correct the distance with the zero Pseudomonas count
Mock.even_cal_Ait_dist[row.names(Mock.even_cal_Ait_dist) == "502-BiomPsy-931", ] <-
  vegdist(t(Mock.even_cal_Ait[!grepl(
    "Pseudo", row.names(Mock.even_cal_Ait)), c("502-BiomPsy-931", "DNA")]),
    method = "euclidean")

Mock.even_cal_Ait_dist <- 
  Mock.even_cal_Ait_dist[row.names(Mock.even_cal_Ait_dist) != "DNA", , drop = F]
colnames(Mock.even_cal_Ait_dist) <- "Ait_cal"
Mock.even_cal_Ait_dist <- merge(
  Mock.even_cal_Ait_dist, mock.even_Ait, by.x = 0, by.y = "Row.names",
  all = TRUE)

# ------------------------------------------------------------------------------
# Statistical test
# ------------------------------------------------------------------------------

Even_Ait_test <- Mock.even_cal_Ait_dist %>%
  rename("Uncorrected" = Aitchison_dist, "Corrected" = Ait_cal) %>%
  select(Dil, Uncorrected, Corrected, Prot) %>%
  melt(.)

# Main, 10^6, p = 0.031
wilcox.test(value ~ variable, data = Even_Ait_test %>% 
              filter(Dil == "10^6 cells", !grepl("Z_T", Prot)), paired = TRUE)
# Main, 10^8, p = 0.031
wilcox.test(value ~ variable, data = Even_Ait_test %>% 
              filter(Dil == "10^8 cells", !grepl("Z_T", Prot)), paired = TRUE)
# Main, 10^5, p = 0.031
wilcox.test(value ~ variable, data = Even_Ait_test %>% 
              filter(Dil == "10^5 cells", !grepl("Z_T", Prot)), paired = TRUE)
# Main, 10^4, p = 0.031
wilcox.test(value ~ variable, data = Even_Ait_test %>% 
              filter(Dil == "10^4 cells", !grepl("Z_T", Prot)), paired = TRUE)

# Supp, 10^6, p = 0.0078
wilcox.test(value ~ variable, data = Even_Ait_test %>% 
              filter(Dil == "10^6 cells"), paired = TRUE)
# Supp, 10^8, p = 0.0078
wilcox.test(value ~ variable, data = Even_Ait_test %>% filter(Dil == "10^8 cells"),
            paired = TRUE)
# Supp, 10^5, p = 0.055
wilcox.test(value ~ variable, data = Even_Ait_test %>% 
              filter(Dil == "10^5 cells"), paired = TRUE)
# Supp, 10^4, p = 0.0078
wilcox.test(value ~ variable, data = Even_Ait_test %>% 
              filter(Dil == "10^4 cells"), paired = TRUE)

################################################################################
#
# Correct staggered mock samples with this bias
#
################################################################################

# Add protocol to observed data
Mock_log_cal <- seqtab_mock.log %>% 
  mutate(DNA = apply(across(S.ID$DNA.log), 1, geom_mean)) %>% 
  select(all_of(S.ID$mock.log), "DNA") %>%
  rownames_to_column("taxon") %>% 
  melt(., variable.name = "Sample_ID") %>%
  mutate(Group = case_when(
    grepl("Pseudo|Salmon|Esche", taxon) ~ "rods, gram neg.",
    grepl("Bacillus|Lactoba", taxon, ignore.case = F) ~ "rods, gram pos.",
    grepl("List|Enteroc|Staph", taxon) ~ "coccoid, gram pos."))

# Merge with metafile (add protocol information)
Mock_log_cal <- merge(Mock_log_cal, Metafile[, c("Sample_ID", "Prot")], 
                           by = "Sample_ID", all.x = TRUE)

# Merge bias (3 groups * 8 prot) with obs data ((8 prots * 2 dilutions+1) * 4 species)
Mock_log_cal <- Mock_log_cal %>%
  merge(., Bias_even_106, by = c("Prot", "Group"), all = TRUE) 

# Calibrate counts, make relative abundances
Mock_log_cal <- Mock_log_cal %>%
  # Calibrate observed counts with bias (divide observed values by bias)
  mutate(val_cal = value/estimate_gm) %>% # estimate_gm!
  # Scale each sample to 100%
  group_by(Sample_ID) %>%
  mutate(val_scale = value/sum(value),
         val_cal_scale = val_cal/sum(val_cal)) %>%
  mutate(BD_BC = case_when(is.na(estimate_gm) ~ val_scale, TRUE ~ val_cal_scale), # scaled data for BC
         BD_Ait = case_when(is.na(estimate_gm) ~ value, TRUE ~ val_cal)) # count data for AIT

# ------------------------------------------------------------------------------
# BC
# ------------------------------------------------------------------------------

# Data in final format: Calibrated DNA mock samples
Mock.log_cal_BC <- Mock_log_cal %>% 
  select(taxon, Sample_ID, BD_BC) %>%
  pivot_wider(names_from = taxon, values_from = BD_BC) %>%
  column_to_rownames("Sample_ID") %>% 
  t(.) %>% as.data.frame()

# Calculate BC distance, select only column "Expected"
Mock.log_cal_BC_dist <- data.frame(BC_cal = as.matrix(
  vegdist(t(Mock.log_cal_BC), method = "bray"))[, "DNA"])
# Remove the distance between DNA and DNA sample
Mock.log_cal_BC_dist <- Mock.log_cal_BC_dist %>% filter(row.names(.) != "DNA")
Mock.log_cal_BC_dist <- merge(
  Mock.log_cal_BC_dist, mock.log_BC, by.x = 0, by.y = "Row.names", all = TRUE) %>%
  rename("Original" = BC_dist, "Corrected" = BC_cal)

# Plot the result - boxplot, 10^8 & 10^6, main figure
annotations <- data.frame(
  x = c(0.5, 0.5),
  y = c(0.038, 0.0045),
  Dil = factor(c("10^8 cells", "10^6 cells"), levels = c("10^8 cells", "10^6 cells")),
  label = c("p = 0.031", "p = 0.031"))
svg("Output/Plots/F06C_Dist-log_2023-08-16.svg", 
       height = 4.2, width = 3)
Mock.log_cal_BC_dist %>% 
  mutate(Prot = factor(Prot, levels = rev(sort(unique(Metafile$Prot))))) %>%
  select(Dil, Original, Corrected, Prot) %>% melt(.) %>%
  filter(!grepl("Z_T", Prot)) %>%
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot() + geom_point() +
  geom_line(aes(group = Prot, colour = Prot), size = 2) +
  scale_y_log10() + ylab("BC distance to DNA mock") +
  facet_grid(Dil ~., scales = "free_y") +
  scale_colour_manual(values = Col_prot) + plot_theme + 
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = rel(1.2), colour = "black")) +
  geom_text(data = annotations, aes(x = x, y = y, label = label), 
            hjust = 0, vjust = 0)
dev.off()

# Plot the result - boxplot, 10^8 & 10^6, supp figure
annotations <- data.frame(
  x = c(0.5, 0.5),
  y = c(0.0275, 0.0045),
  Dil = factor(c("10^8 cells", "10^6 cells"), levels = c("10^8 cells", "10^6 cells")),
  label = c("p = 0.20", "p = 0.20"))
svg("Output/Plots/S07C_Dist-log_2023-08-16.svg", 
       height = 4.2, width = 3)
Mock.log_cal_BC_dist %>% 
  # To have yellow and orange line in the background
  mutate(Prot = factor(Prot, levels = rev(sort(unique(Metafile$Prot))))) %>%
  select(Dil, Original, Corrected, Prot) %>% melt(.) %>%
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot() + geom_point() +
  geom_line(aes(group = Prot, colour = Prot), size = 2) +
  scale_y_log10() + ylab("BC distance to DNA mock") +
  facet_grid(Dil ~., scales = "free_y") +
  scale_colour_manual(values = Col_prot) + plot_theme + 
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = rel(1.2), colour = "black")) +
  geom_text(data = annotations, aes(x = x, y = y, label = label), 
            hjust = 0, vjust = 0)
dev.off()

# ------------------------------------------------------------------------------
# Statistical test - Bray-Curtis
# ------------------------------------------------------------------------------

Log_BC_test <- Mock.log_cal_BC_dist %>%
  select(Dil, Original, Corrected, Prot) %>% melt(.)

# Main, 10^8, p = 0.031
wilcox.test(value ~ variable, data = filter(Log_BC_test, Dil == "10^8 cells",
                                            !grepl("Z_T", Prot)), paired = TRUE)
# Main, 10^6, p = 0.031
wilcox.test(value ~ variable, data = filter(Log_BC_test, Dil == "10^6 cells",
                                            !grepl("Z_T", Prot)), paired = TRUE)

# Supp, 10^8, p = 0.20
wilcox.test(value ~ variable, data = Log_BC_test %>% filter(Dil == "10^8 cells"),
            paired = TRUE)
# Supp, 10^6, p = 0.20
wilcox.test(value ~ variable, data = Log_BC_test %>% filter(Dil == "10^6 cells"),
            paired = TRUE)

# ------------------------------------------------------------------------------
# Aitchison
# ------------------------------------------------------------------------------

# Calculate Aitchison distance
Mock.log_cal_Ait <- Mock_log_cal %>%
  select(Sample_ID, taxon, BD_Ait) %>%
  pivot_wider(names_from = taxon, values_from = BD_Ait) %>%
  column_to_rownames("Sample_ID") %>% 
  t(.) %>% as.data.frame()
# Aitchison distance
Mock.log_cal_Ait <- apply(Mock.log_cal_Ait, 2, function(x) clr(x))
round(colSums(Mock.log_cal_Ait), 10) # check
# Calculate Aitchison distance, select only column "Expected"
Mock.log_cal_Ait <- data.frame(Aitchison_dist = as.matrix(
  vegdist(t(Mock.log_cal_Ait), method = "euclidean"))[, "DNA"])

Mock.log_cal_Ait <- 
  Mock.log_cal_Ait[row.names(Mock.log_cal_Ait) != "DNA", , drop = F]
colnames(Mock.log_cal_Ait) <- "Ait_cal"
Mock.log_cal_Ait <- merge(
  Mock.log_cal_Ait, mock.log_Ait, by.x = 0, by.y = "Row.names",
  all = TRUE)

# ------------------------------------------------------------------------------
# Statistical test
# ------------------------------------------------------------------------------

Log_Ait_test <- Mock.log_cal_Ait %>%
  rename("Uncorrected" = Aitchison_dist, "Corrected" = Ait_cal) %>%
  select(Dil, Uncorrected, Corrected, Prot) %>%
  melt(.)

# Main, 10^8, p = 0.031
wilcox.test(value ~ variable, data = filter(Log_Ait_test, Dil == "10^8 cells",
                                            !grepl("Z_T", Prot)), paired = TRUE)
# Main, 10^6, p = 0.031
wilcox.test(value ~ variable, data = filter(Log_Ait_test, Dil == "10^6 cells",
                                            !grepl("Z_T", Prot)), paired = TRUE)

# Supp, 10^8, p = 0.31
wilcox.test(value ~ variable, data = Log_Ait_test %>% filter(Dil == "10^8 cells"),
            paired = TRUE)
# Supp, 10^6, p = 0.64
wilcox.test(value ~ variable, data = Log_Ait_test %>% filter(Dil == "10^6 cells"),
            paired = TRUE)

################################################################################
#
# Correct biases in pure spike samples
#
################################################################################

# Define observed data, only cell samples (no DNA samples)
Obs_spike_plus <- Obs_spike
Obs_spike_plus[Obs_spike_plus == 0] <- 0.5

# Add protocol to observed data
Obs_spike_plus <- Obs_spike_plus %>% 
  # Add expected data
  merge(., Expected_spike, by.x = 0, by.y = "Species", all = TRUE) %>% 
  rename(taxon = Row.names, Expected = Abund_16S)

# Calculate BC distance between observed and DNA mock
Obs_spike_plus_BC <- Obs_spike_plus %>% column_to_rownames("taxon") %>%
  mutate(across(everything(), ~ .x/sum(.x)))
Obs_spike_plus_BC <- data.frame(BC_dist = as.matrix(
  vegdist(t(Obs_spike_plus_BC), method = "bray"))[, "Expected"]) %>%
  filter(row.names(.) != "Expected")

# Calculate Aitchison distance between observed and DNA mock
Obs_spike_plus_Ait <- Obs_spike_plus %>% column_to_rownames("taxon")
# Aitchison distance
Obs_spike_plus_Ait <- apply(Obs_spike_plus_Ait, 2, function(x) clr(x))
round(colSums(Obs_spike_plus_Ait), 10) # check
Obs_spike_plus_Ait <- data.frame(Aitchison_dist = as.matrix(
  vegdist(t(Obs_spike_plus_Ait), method = "euclidean"))[, "Expected"]) %>%
  filter(row.names(.) != "Expected")

# Continue as usual 
Obs_spike_plus <- Obs_spike_plus %>%
  melt(., variable.name = "Sample_ID") %>% 
  mutate(Group = case_when(grepl("Imt", taxon) ~ "rods, gram neg.",
                           grepl("Allo", taxon) ~ "rods, gram pos.",
                           grepl("True", taxon) ~ "coccoid, gram pos.")) %>%
  # Merge with metafile
  merge(., Metafile[, c("Sample_ID", "Prot")], by = "Sample_ID", all.x = TRUE)

# Merge bias (3 groups *8 prot) with obs data (3 species * (15+1) samples)
Bias_spike_cal <- Bias_even_106 %>%
  merge(., Obs_spike_plus, by = c("Prot", "Group"), all = TRUE) %>% 
  # Calibrate observed counts with bias (divide observed values by bias)
  mutate(val_cal = value/estimate_gm) %>% 
  # Scale each sample to 100%
  group_by(Sample_ID) %>%
  mutate(val_scale = value/sum(value),
         val_cal_scale = val_cal/sum(val_cal)) %>%
  mutate(BD_BC = case_when(is.na(estimate_gm) ~ val_scale, TRUE ~ val_cal_scale), # scaled data for BC
         BD_Ait = case_when(is.na(estimate_gm) ~ value, TRUE ~ val_cal)) # count data for AIT

# ------------------------------------------------------------------------------
# BC
# ------------------------------------------------------------------------------

# Data in final format: Calibrated mock samples (calibrierte rel. ab. in wide format)
Bias_spike_cal_BC <- Bias_spike_cal %>% 
  select(taxon, Sample_ID, BD_BC) %>%
  pivot_wider(names_from = taxon, values_from = BD_BC) %>%
  column_to_rownames("Sample_ID") %>% 
  t(.) %>% as.data.frame()

# Calculate BC distance, select only column "Expected"
BC_spike_cal <- data.frame(BC_cal = as.matrix(
  vegdist(t(Bias_spike_cal_BC), method = "bray"))[, "Expected"])
# Correct the zeros

# Remove the distance between Expected and Expected sample
BC_spike_cal <- 
  BC_spike_cal[!row.names(BC_spike_cal) %in% c("Expected"), , drop = F]

# Merge BC distances (original & calibrated), and metafile
BC_spike_cal <- BC_spike_cal %>%
  merge(., Obs_spike_plus_BC, by = 0, all = TRUE) %>% 
  merge(Metafile, ., by.x = "Sample_ID", by.y = "Row.names", all.y = TRUE) %>%
  rename("Original" = BC_dist, "Corrected" = BC_cal)

# Plot the result - boxplot, main figure
annotations <- data.frame(
  x = c(0.5, 0.5),
  y = c(0.05, 0.008),
  Dil = c("10^5 cells", "6*10^3 cells"),
  label = c("p = 0.063", "p = 0.031"))
svg("Output/Plots/F06D_Dist-spike_2023-08-16.svg", 
       height = 4.2, width = 3)
BC_spike_cal %>% 
  # Important for line order
  mutate(Prot = factor(Prot, levels = rev(sort(unique(Metafile$Prot))))) %>%
  select(Dil, Original, Corrected, Prot) %>% melt(.) %>%
  filter(!grepl("Z_T", Prot)) %>%
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot() + geom_point() +
  geom_line(aes(group = Prot, colour = Prot), size = 2) +
  scale_y_log10() + ylab("BC distance to expected composition") +
  facet_grid(Dil ~., scales = "free_y") +
  scale_colour_manual(values = Col_prot) + plot_theme + 
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = rel(1.2), colour = "black")) +
  geom_text(data = annotations, aes(x = x, y = y, label = label), 
            hjust = 0, vjust = 0)
dev.off()

# Plot the result - boxplot, supp figure
annotations <- data.frame(
  x = c(0.5, 0.5),
  y = c(0.05, 0.008),
  Dil = c("10^5 cells", "6*10^3 cells"),
  label = c("p = 0.039", "p = 0.016"))
svg("Output/Plots/S07D_Dist-spike_2023-08-16.svg", 
       height = 4.2, width = 3)
BC_spike_cal %>% 
  # Important for line order
  mutate(Prot = factor(Prot, levels = rev(sort(unique(Metafile$Prot))))) %>%
  select(Dil, Original, Corrected, Prot) %>% melt(.) %>% 
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot() + geom_point() +
  geom_line(aes(group = Prot, colour = Prot), size = 2) +
  scale_y_log10() + ylab("BC distance to expected composition") +
  facet_grid(Dil ~., scales = "free_y") +
  scale_colour_manual(values = Col_prot) + plot_theme + 
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = rel(1.2), colour = "black")) +
  geom_text(data = annotations, aes(x = x, y = y, label = label), 
            hjust = 0, vjust = 0)
dev.off()

# ------------------------------------------------------------------------------
# Statistical test - Bray-Curtis
# ------------------------------------------------------------------------------

BC_test3 <- BC_spike_cal %>%
  select(Dil, Original, Corrected, Prot) %>% melt(.)

# Main, 10^5, p = 0.063
wilcox.test(value ~ variable, data = filter(BC_test3, Dil == "10^5 cells",
                                            !grepl("Z_T", Prot)), paired = TRUE)
# Main, 6*10^3, p = 0.031
wilcox.test(value ~ variable, data = filter(BC_test3, Dil == "6*10^3 cells",
                                            !grepl("Z_T", Prot)), paired = TRUE)

# Supp, 10^5, p = 0.039
wilcox.test(value ~ variable, data = BC_test3 %>% filter(Dil == "10^5 cells"),
            paired = TRUE)
# Supp, 6*10^3, p = 0.016
wilcox.test(value ~ variable, data = filter(BC_test3, Dil == "6*10^3 cells"), 
            paired = TRUE)

# ------------------------------------------------------------------------------
# Aitchison
# ------------------------------------------------------------------------------

# Data in final format: Calibrated mock samples (calibrierte rel. ab. in wide format)
Aitchison_spike_cal <- Bias_spike_cal %>% 
  select(taxon, Sample_ID, BD_Ait) %>%
  pivot_wider(names_from = taxon, values_from = BD_Ait) %>%
  column_to_rownames("Sample_ID") %>% 
  t(.) %>% as.data.frame()
# Aitchison distance
Aitchison_spike_cal <- apply(Aitchison_spike_cal, 2, function(x) clr(x))
round(colSums(Aitchison_spike_cal), 10) # check
# Calculate Aitchison distance, select only column "Expected"
Aitchison_spike_cal <- data.frame(Ait_cal = as.matrix(
  vegdist(t(Aitchison_spike_cal), method = "euclidean"))[, "Expected"]) %>%
  filter(row.names(.) != "Expected")

# Merge BC distances (original & calibrated), and metafile
Aitchison_spike_cal <- Aitchison_spike_cal %>%
  merge(., Obs_spike_plus_Ait, by = 0, all = TRUE) %>% 
  merge(Metafile, ., by.x = "Sample_ID", by.y = "Row.names", all.y = TRUE) %>%
  rename("Original" = Aitchison_dist, "Corrected" = Ait_cal)

# ------------------------------------------------------------------------------
# Statistical test
# ------------------------------------------------------------------------------

Ait_test3 <- Aitchison_spike_cal %>%
  select(Dil, Original, Corrected, Prot) %>%
  melt(.)

# Main, 10^5, p = 0.031
wilcox.test(value ~ variable, data = filter(Ait_test3, Dil == "10^5 cells",
                                            !grepl("Z_T", Prot)), paired = TRUE)
# Supp, 10^5, p = 0.016
wilcox.test(value ~ variable, data = Ait_test3 %>% filter(Dil == "10^5 cells"),
            paired = TRUE)
# Main, 6*10^3, p = 0.063
wilcox.test(value ~ variable, data = filter(Ait_test3, Dil == "6*10^3 cells",
                                            !grepl("Z_T", Prot)), paired = TRUE)
# Supp, 6*10^3, p = 0.078
wilcox.test(value ~ variable, data = Ait_test3 %>% filter(Dil == "6*10^3 cells"),
            paired = TRUE)

# ------------------------------------------------------------------------------
# Plot legend
# ------------------------------------------------------------------------------

# Legend - main figure
svg("Output/Plots/F06Legend_2023-08-16.svg", height = 2.5, width = 1)
p <- BC_spike_cal %>% 
  mutate(Prot = factor(Prot, levels = rev(sort(unique(Metafile$Prot))))) %>%
  filter(!grepl("Z_T", Prot)) %>%
  select(Dil, Original, Corrected, Prot) %>%
  melt(.) %>%
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot() + geom_point() +
  geom_line(aes(group = Prot, colour = Prot), size = 2) +
  scale_y_log10() + ylab("BC distance to DNA mock") +
  facet_grid(Dil ~., scales = "free_y") +
  scale_colour_manual("Protocol", values = rev(Col_prot[1:8])) +
  plot_theme + 
  guides(colour = guide_legend(reverse = TRUE)) +
  theme(legend.position = "left",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = rel(1.2), colour = "black"))
plot_grid(as_ggplot(get_legend(p)))
dev.off()

# Legend - supp figure
svg("Output/Plots/S07Legend_2023-08-16.svg", height = 2.5, width = 1)
p <- BC_spike_cal %>% 
  mutate(Prot = factor(Prot, levels = rev(sort(unique(Metafile$Prot))))) %>%
  select(Dil, Original, Corrected, Prot) %>%
  melt(.) %>%
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot() + geom_point() +
  geom_line(aes(group = Prot, colour = Prot), size = 2) +
  scale_y_log10() + ylab("BC distance to DNA mock") +
  facet_grid(Dil ~., scales = "free_y") +
  scale_colour_manual("Protocol", values = rev(Col_prot[1:8])) +
  plot_theme + 
  guides(colour = guide_legend(reverse = TRUE)) +
  theme(legend.position = "left",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = rel(1.2), colour = "black"))
plot_grid(as_ggplot(get_legend(p)))
dev.off()

# ------------------------------------------------------------------------------
# Median bias reduction (without Z_T)
# ------------------------------------------------------------------------------

# Even mock (from 0.18 to 0.097 BC, decrease by a factor of 1.75 or 43%)
Mock.even_cal_BC_dist %>%
  filter(!grepl("Z_T", Prot)) %>%
  mutate(Fac = Original/Corrected,
         Perc_dec = (Original-Corrected)/Original*100,
         Med_perc = median(Perc_dec),
         Med_fac = median(Fac),
         Med_ori = median(Original),
         Med_cor = median(Corrected))

# Staggered mock (from 0.14 to 0.037, decrease by a factor of 3.76 or 73%)
Mock.log_cal_BC_dist %>%
  filter(!grepl("Z_T", Prot)) %>%
  mutate(Fac = Original/Corrected,
         Perc_dec = (Original-Corrected)/Original*100,
         Med_perc = median(Perc_dec),
         Med_fac = median(Fac),
         Med_ori = median(Original),
         Med_cor = median(Corrected))

# Spike-in mock (from 0.29 to 0.11, decrease by a factor of 2.36 or 57%)
BC_spike_cal %>% 
  filter(!grepl("Z_T", Prot)) %>%
  mutate(Fac = Original/Corrected,
         Perc_dec = (Original-Corrected)/Original*100,
         Med_perc = median(Perc_dec),
         Med_fac = median(Fac),
         Med_ori = median(Original),
         Med_cor = median(Corrected))

# Even & staggered together (from 0.18 to 0.08 BC, decrease by a factor of 2.12 or 53%)
Mock.even_cal_BC_dist %>%
  rbind.data.frame(., Mock.log_cal_BC_dist) %>% 
  rename(Sample_ID = Row.names) %>%
  #rbind.data.frame(., BC_spike_cal) %>% 
  filter(!grepl("Z_T", Prot)) %>%
  mutate(Fac = Original/Corrected,
         Perc_dec = (Original-Corrected)/Original*100,
         Med_perc = median(Perc_dec),
         Med_fac = median(Fac),
         Med_ori = median(Original),
         Med_cor = median(Corrected))

# BC distance between DNA mocks and expected composition - EVEN mock
seqtab_mock.even[, S.ID$DNA.even] %>%
  merge(., Expected_mock[, c("Species", "Abund_even_16S")], by.x = 0, 
        by.y = "Species") %>%
  mutate(across(where(is.numeric), ~ .x/sum(.x))) %>%
  column_to_rownames("Row.names") %>% 
  t(.) %>%
  vegdist(., method = "bray") %>% as.matrix()
# BC distance between DNA mocks and expected composition - STAGGERED mock
seqtab_mock.log[, S.ID$DNA.log] %>%
  merge(., Expected_mock[, c("Species", "Abund_log_16S")], by.x = 0, 
        by.y = "Species") %>%
  mutate(across(where(is.numeric), ~ .x/sum(.x))) %>%
  column_to_rownames("Row.names") %>% 
  t(.) %>%
  vegdist(., method = "bray") %>% as.matrix()  
median(c(0.06826205, 0.05474897, 0.044155124, 0.040652443)) # 0.049

# ------------------------------------------------------------------------------
# Median distance between original dilutions
# ------------------------------------------------------------------------------

# BC distance between DNA mocks and expected composition - EVEN mock
dist_betw_dil <- c()
for (i in c(13:20)) { # check: S.ID[c(13:36)]
  dist_betw_dil <- append(
    dist_betw_dil, 
    seqtab_mock.even[, S.ID[[i]]] %>%
      mutate(across(where(is.numeric), ~ .x/sum(.x))) %>% t(.) %>%
      vegdist(., method = "bray", upper = FALSE) %>% as.vector())
}
# STAGGERED mock
for (i in c(21:28)) { # check: S.ID[c(13:36)]
  dist_betw_dil <- append(
    dist_betw_dil, 
    seqtab_mock.log[, S.ID[[i]]] %>%
      mutate(across(where(is.numeric), ~ .x/sum(.x))) %>% t(.) %>%
      vegdist(., method = "bray", upper = FALSE) %>% as.vector())
}
median(dist_betw_dil)
