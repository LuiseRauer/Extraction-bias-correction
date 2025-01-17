################################################################################
#
# Zymo-Qiagen comparison mock analysis
# Detailed analysis of sequence errors, chimeras, and unknowns
#
################################################################################

# Load required packages
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(scales)

###save(list = ls(), file = "R_QZmock_scripts1-6_2024-08-28.RData")
###load("R_QZmock_scripts1-6_2024-08-28.RData", verbose = TRUE)

# Calculate relative abundances from ASV counts
CCC_detail <- CCC_res %>%
  mutate(across(S.ID$all, ~ .x/sum(.x)))

################################################################################
#
# Sequence errors
#
################################################################################

# Create an ASV table with only sequence errors 
CCC_detail_err <- CCC_detail %>%
  # Keep only sequence errors
  filter(Category == "Amplification_error") %>% 
  # Select all mock samples
  select(ASV_ID, Min_LV_dist, LV_8, all_of(c(S.ID$mock.all))) %>%
  # Tranform to long format
  melt(., id.vars = c("ASV_ID", "Min_LV_dist", "LV_8"), 
       variable.name = "Sample_ID", value.name = "Rel_abund") %>%
  # Merge with Metafile
  merge(., Metafile[, c("Sample_ID", "Prot", "Broad_type", "Dil")], 
        by = "Sample_ID") %>%
  # Remove zero relative abundances for better computing time
  filter(Rel_abund != 0)

# Determine ASV order in the plot by decreasing span and rel. abundance
ASV_order <- CCC_detail_err %>%
  group_by(Broad_type, Dil, ASV_ID) %>%
  # Calculate ASV span per dilution and mock type (min 1, max 8)
  mutate(Span = sum(Rel_abund > 0)) %>%
  ungroup() %>%
  # Order by descending span and descending relative abundance
  arrange(Broad_type, Dil, desc(Span), desc(Rel_abund)) %>% 
  select(ASV_ID) %>% unlist(.) %>% unique(.)

# Add data for x axis formatting and order samples and ASVs
CCC_detail_err <- CCC_detail_err %>% 
  # Add NA data for samples without any sequence errors for visualization
  add_row(ASV_ID = "ASV_0106", LV_8 = "Bacillus_subtilis", Rel_abund = NA,
          Min_LV_dist = 8, 
          Prot = c("Q_S_z", "Z_T_q", "Z_T_q", "Q_T_z", "Z_T_q"), 
          Dil = paste0(c("10^6", "10^6", "10^5", "6*10^3", "6*10^3"), " cells"), 
          Broad_type = c("Staggered mock", "Staggered mock",
                         "Spike-in", "Spike-in", "Spike-in")) %>%
  # Order ASVs, species facets, and samples in mocks and dilutions
  mutate(ASV_ID = factor(ASV_ID, levels = rev(ASV_order)),
         LV_8 = factor(gsub("_", " ", LV_8), levels = gsub(
           "_", " ", unique(Table_16S_279$ASV.ID))),
         Broad_type = factor(Broad_type, levels = c(
           "Even mock", "Staggered mock", "Spike-in")),
         Dil = factor(Dil, levels = c(
           paste0("10^", c(8, 6, 5, 4), " cells"), "6*10^3 cells", "D")))

# Plot the data
svg("Output/Plots/S01_CCC-sequence-errors_2023-05-26.svg",
       width = 7.5, height = 7)
ggplot(CCC_detail_err, aes(x = Prot, y = ASV_ID, size = Rel_abund, 
                           colour = as.character(Min_LV_dist))) +
  geom_point() +
  ggh4x::facet_nested(LV_8 ~ Broad_type + Dil, scales = "free", space = "free") +
  scale_size("Relative abundance", range = c(0.3, 2),
             breaks = c(0.04, 0.08, 0.12, 0.16)) +
  scale_colour_manual("LV distance", values = c(Col_class[2:9])) +
  xlab("Protocol") + ylab("ASV") + plot_theme +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, 
                                    face = "italic", size = rel(0.7)),
        strip.text.x = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(0.5)),
        axis.text.x = element_text(angle = 90, size = rel(0.5),
                                   hjust = 0, vjust = 0.5),
        panel.spacing = unit(0, "lines"))
dev.off()

# Text: proportion of Lact/Esch/Salm reads of all sequence errors in even mock
CCC_detail_err %>% 
  # select cell and DNA samples of the even mock 
  filter(Sample_ID %in% c(S.ID$mock.even, S.ID$DNA.even)) %>%
  mutate(Err_group = case_when(grepl("Lacto|Esche|Salm", LV_8) ~ "Major",
                               TRUE ~ "Minor")) %>%
  # Group by species' category and sample ID
  group_by(Err_group, Sample_ID) %>% 
  # Calculate rel. abundance of the three taxa of all sequence errors per sample
  summarise(val_sum = sum(Rel_abund, na.rm = T)) %>%
  # Transform the data into "wide" format
  spread(Err_group, val_sum) %>% 
  # Replace NA by 0 for calculating the median
  replace_na(list(Minor = 0)) %>% 
  mutate(Prop = Major/Minor) %>% summarise(med = median(Prop)) # 76.1 %

################################################################################
#
# Chimeras
#
################################################################################

# Are there chimeras in samples with 10^5 and 10^4 cells and spike-in?
CCC_detail %>%
  filter(Category == "Mock_chimera") %>% 
  mutate(Sum_abund_lowinput = rowSums(across(all_of(c(S.ID$spike, S.ID$m.e.dil3, 
                                                      S.ID$m.e.dil4))))) %>%
  filter(Sum_abund_lowinput > 0) %>% 
  select(ASV_ID, Sum_abund_lowinput) # Only 4 chimeras, max. rel. abund. 0.01% 

# ------------------------------------------------------------------------------
# Which species are forming chimeras in 10^8 and 10^6 samples?
# ------------------------------------------------------------------------------

# Detailed analysis of chimeras in 10^8, 10^6 and DNA mock samples
CCC_detail_chim <- CCC_detail %>%
  filter(Category == "Mock_chimera") %>% 
  select(ASV_ID, Min_LV_dist, Taxon_L1, Taxon_L2, Taxon_L3, 
         all_of(c(S.ID$m.e.dil1, S.ID$m.e.dil2, S.ID$m.l.dil1, S.ID$m.l.dil2, 
                  S.ID$DNA.even, S.ID$DNA.log))) %>%
  # Shorten taxon names
  mutate(
    tax1 = paste0(substr(Taxon_L1, 1, 1), ". ", sub(".*_", "", Taxon_L1)),
    tax2 = paste0(substr(Taxon_L2, 1, 1), ". ", sub(".*_", "", Taxon_L2)),
    tax3 = case_when(
      is.na(Taxon_L3) ~ NA_character_,
      TRUE ~ paste0(substr(Taxon_L3, 1, 1), ". ", sub(".*_", "", Taxon_L3))))

# Check that no chimera is composed of 2 or 3 identical species
CCC_detail_chim %>%
  mutate(Tax_sum = case_when(Taxon_L1 == Taxon_L2 & is.na(Taxon_L3) ~ "eq",
                             Taxon_L1 == Taxon_L2 & Taxon_L2 == Taxon_L3 ~ "eq", 
                             TRUE ~ "diff")) %>% 
  filter(Tax_sum == "eq")

# Paste shortened and alphabetically ordered taxon names into one column
CCC_detail_chim["Chim_comb"] <- 
  apply(CCC_detail_chim[, c("tax1", "tax2", "tax3")], 1, function(x) 
    paste(sort(x), collapse = ", "))

# Format the data for plotting
CCC_detail_chim <- CCC_detail_chim %>%
  # Keep only the pasted taxon column
  select(-c(Taxon_L1, Taxon_L2, Taxon_L3, tax1, tax2, tax3)) %>%
  # Transform to long format
  melt(., id.vars = c("ASV_ID", "Min_LV_dist", "Chim_comb"), 
       variable.name = "Sample_ID", value.name = "Rel_abund") %>%
  filter(Rel_abund != 0) %>% 
  group_by(Chim_comb, Sample_ID) %>% 
  # Sum relative abundance and count ASVs per sample
  summarise(Rel_abund = sum(Rel_abund),
            Count = n()) %>%
  # Merge with Metafile
  merge(., Metafile[, c("Sample_ID", "Prot", "Broad_type", "Dil")], 
        by = "Sample_ID") %>%
  # Classify into bimeras and trimeras
  mutate(Type = case_when(
    str_count(Chim_comb, "\\. ") == 2 ~ "Bimera", TRUE ~ "Trimera")) %>%
  # Add NA data for protocols without any chimeras
  add_row(Chim_comb = "B. subtilis, E. coli", Type = "Bimera", Dil = "10^6 cells",
          Prot = c("Q_S_z", "Q_T_q", "Q_S_q", "Q_T_q", "Q_T_z"),
          Broad_type = c(rep("Even mock", 2), rep("Staggered mock", 3)),
          Count = 1, Rel_abund = NA) %>%
  # Order dilutions
  mutate(Dil = factor(Dil, levels = c("10^8 cells", "10^6 cells", "D")))

# Plot the result
svg("Output/Plots/S02_CCC-chimeras_2023-04-28.svg",
       width = 7.5, height = 6)
CCC_detail_chim %>%
  mutate(Chim_comb = gsub("_", " ", Chim_comb)) %>%
  ggplot(., aes(x = Prot, y = Chim_comb, size = Rel_abund, 
                colour = as.character(Count))) +
  geom_point() +
  scale_colour_manual("Number of ASVs", values = rev(Col_ASV)) +
  scale_y_discrete(limits = rev) +
  scale_size("Relative abundance", range = c(0.3, 2.5)) +
  ggh4x::facet_nested(Type ~ Broad_type + Dil, scales = "free", space = "free") +
  xlab("Protocol") + 
  ylab("Taxon combination (alphabetically ordered)") + 
  plot_theme +
  theme(axis.text.y = element_text(face = "italic", size = rel(0.7)),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, 
                                   size = rel(0.7)),
        panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = rel(0.7)))
dev.off()

# ------------------------------------------------------------------------------
# Sequence distance and chimera occurrence in 10^8 even and staggered samples
# ------------------------------------------------------------------------------

# Order species by decreasing abundance in the staggered mock
tax_order <- Expected_mock$Species[order(Expected_mock$Abund_log_16S, 
                                         decreasing = TRUE)]
tax_order2 <- sapply(tax_order, function(x) 
  gsub(".*_", paste0(substring(x, 1, 1), ". "), x))

# Select chimeras in the even mock
CCC_detail_chim_even <- CCC_detail %>%
  filter(Category == "Mock_chimera") %>% 
  select(ASV_ID, Min_LV_dist, Taxon_L1, Taxon_L2, Taxon_L3, 
         all_of(S.ID$m.e.dil1))
# Sort the two taxa by decreasing abundance in the staggered mock (first taxon)
tax1 <- apply(CCC_detail_chim_even, 1, function(x) {
  # If 1st and 2nd taxon are identical, order 2nd and 3rd taxon by vector 
  case_when(x["Taxon_L1"] == x["Taxon_L2"] ~ 
              tax_order[min(which(tax_order == x["Taxon_L2"]), 
                            which(tax_order == x["Taxon_L3"]))],
            # else order 1st and 2nd taxon by vector
            TRUE ~ 
              tax_order[min(which(tax_order == x["Taxon_L1"]), 
                            which(tax_order == x["Taxon_L2"]))])
})
# Sort the two taxa by decreasing abundance in the staggered mock (second taxon)
tax2 <- apply(CCC_detail_chim_even, 1, function(x) {
  # If 1st and 2nd taxon are identical, order 2nd and 3rd taxon by vector 
  case_when(x["Taxon_L1"] == x["Taxon_L2"] ~
              tax_order[max(which(tax_order == x["Taxon_L2"]), 
                            which(tax_order == x["Taxon_L3"]))],
            # else order 1st and 2nd taxon by vector
            TRUE ~ 
              tax_order[max(which(tax_order == x["Taxon_L1"]), 
                            which(tax_order == x["Taxon_L2"]))])
}) 
CCC_detail_chim_even <- cbind.data.frame(CCC_detail_chim_even, tax1, tax2) %>%
  # Summarise identical chimera annotations
  group_by(tax1, tax2) %>%
  summarise(across(S.ID$m.e.dil1, sum)) %>%
  ungroup(.) %>%
  # Calculate mean rel. abundance and span per chimera type over protocols
  mutate(mean_ab = rowMeans(across(all_of(S.ID$m.e.dil1))),
         span = rowSpan(across(all_of(S.ID$m.e.dil1)))) %>%
  filter(mean_ab > 0) %>%
  select(tax1, tax2, mean_ab, span) %>%
  mutate(tax1 = factor(tax1, levels = unique(Table_16S_279$ASV.ID)[1:8]),
         tax2 = factor(tax2, levels = unique(Table_16S_279$ASV.ID)[1:8])) %>% 
  ungroup() %>%
  # Add empty taxon combinations
  complete(tax1, tax2) 

# Staggered mock
CCC_detail_chim_log <- CCC_detail %>%
  filter(Category == "Mock_chimera") %>% 
  select(ASV_ID, Min_LV_dist, Taxon_L1, Taxon_L2, Taxon_L3, 
         all_of(S.ID$m.l.dil1))  
# Sort the two taxa anti-alphabetically (switch order of tax1 and tax2)
tax2 <- apply(CCC_detail_chim_log, 1, function(x) {
  case_when(x["Taxon_L1"] == x["Taxon_L2"] ~ 
              tax_order[min(which(tax_order == x["Taxon_L2"]), 
                            which(tax_order == x["Taxon_L3"]))],
            TRUE ~ 
              tax_order[min(which(tax_order == x["Taxon_L1"]), 
                            which(tax_order == x["Taxon_L2"]))])
})
tax1 <- apply(CCC_detail_chim_log, 1, function(x) {
  case_when(x["Taxon_L1"] == x["Taxon_L2"] ~
              tax_order[max(which(tax_order == x["Taxon_L2"]), 
                            which(tax_order == x["Taxon_L3"]))],
            TRUE ~ 
              tax_order[max(which(tax_order == x["Taxon_L1"]), 
                            which(tax_order == x["Taxon_L2"]))])
})
CCC_detail_chim_log <- cbind.data.frame(CCC_detail_chim_log, tax1, tax2) %>%
  # Summarise identical chimera annotations
  group_by(tax1, tax2) %>%
  summarise(across(S.ID$m.l.dil1, sum)) %>%
  ungroup(.) %>%
  # Calculate mean rel. abundance and span per chimera type over protocols
  mutate(mean_ab = rowMeans(across(all_of(S.ID$m.l.dil1))),
         span = rowSpan(across(all_of(S.ID$m.l.dil1)))) %>% 
  filter(mean_ab > 0) %>%
  select(tax1, tax2, mean_ab, span) %>%
  mutate(tax1 = factor(tax1, levels = unique(Table_16S_279$ASV.ID)[1:8]),
         tax2 = factor(tax2, levels = unique(Table_16S_279$ASV.ID)[1:8])) %>% 
  ungroup() %>%
  # Add empty taxon combinations
  complete(tax1, tax2) 

# Combine results from even and staggered mock
CCC_detail_chim_plot <- CCC_detail_chim_even %>% 
  mutate(Mock = "even") %>%
  rbind.data.frame(., CCC_detail_chim_log %>% mutate(Mock = "log")) %>% 
  replace_na(list(span = 1)) %>%
  mutate(tax1 = factor(gsub("_", " ", tax1), levels = gsub("_", " ", tax_order)),
         tax2 = factor(gsub("_", " ", tax2), levels = gsub("_", " ", tax_order))
  ) %>% 
  merge(melt(Exp_dist_279_mean, value.name = "LV"), by.x = c("tax1", "tax2"), 
        by.y = c("Var1", "Var2"), all.x = TRUE)

svg("Output/Plots/F03B_CCC-chimeras_2023-08-06.svg", width = 5.09, height = 3.26)
CCC_detail_chim_plot %>%
  ggplot(., aes(x = as.numeric(tax1), y = as.numeric(tax2), 
                size = as.numeric(mean_ab), colour = LV)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  scale_colour_gradientn(
    "LV sequence distance", limits = c(0, 150), 
    breaks = seq(0, 140, 35),
    colours = c("#1E47C2", "#226BCA", "#2288BF", 
                "#17BF71", "#B0E03F", "#FFED47",
                "#F5C149", "#EB954B", "#E1694D", "#D63D4E", "#B12F3E")) +
  scale_size("Chimera rel. abund.", breaks = c(0.0002, 0.002, 0.02), 
             labels = c("0.0002", "0.002", "0.02"), range = c(1, 5)) +
  scale_x_continuous(breaks = seq(1, 8, 1), labels = tax_order2,
                     sec.axis = dup_axis()) +
  scale_y_continuous(breaks = seq(1, 8, 1), labels = tax_order2,
                     sec.axis = dup_axis()) +
  plot_theme +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                          face = "italic"),
        axis.text.x.top = element_blank(),
        axis.text.y.left = element_text(face = "italic"),
        axis.text.y.right = element_blank(),
        axis.title = element_text(size = rel(0.9)),
        plot.title.position = "plot",
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank(),
        plot.margin = margin(c(0, 150, 0, 0)),
        legend.position = c(1.65, 0.18),
        axis.ticks = element_blank()) +
  ggtitle("                Chimera formation by taxon") +
  xlab("Even mock") + ylab("Staggered mock") +
  guides(size = guide_legend(order = 2, reverse = TRUE), 
         colour = guide_colorbar(order = 1))
dev.off()

# ------------------------------------------------------------------------------
# Correlation of chimera proportion & the number of input cells per sample
# ------------------------------------------------------------------------------

# Calculate the proportion of chimeras per mock sample
CCC_detail_chim_corr <- CCC_detail %>%
  filter(Category == "Mock_chimera") %>% 
  select(all_of(S.ID$mock.all)) %>%
  colSums(.) %>% as.data.frame(.) %>%
  merge(Metafile, by.x = 0, by.y = "Sample_ID")
# Replace dilution values for plotting
CCC_detail_chim_corr$Dil <- as.character(CCC_detail_chim_corr$Dil)
CCC_detail_chim_corr$Dil[CCC_detail_chim_corr$Dil == "D"] <- "1.3e7"
CCC_detail_chim_corr$Dil <- gsub(" cells", "", CCC_detail_chim_corr$Dil)
CCC_detail_chim_corr$Dil <- gsub("0\\^", "e", CCC_detail_chim_corr$Dil)
CCC_detail_chim_corr$Dil <- gsub("6\\*1", "6", CCC_detail_chim_corr$Dil)
CCC_detail_chim_corr$Dil <- as.numeric(CCC_detail_chim_corr$Dil)
colnames(CCC_detail_chim_corr)[colnames(CCC_detail_chim_corr) == "."] <- "Chim_prop"

# Spearman correlation between input cell number and chimera read proportion
cor.test(CCC_detail_chim_corr$Chim_prop, CCC_detail_chim_corr$Dil, 
         method = "spearman")

# Text: maximum proportion of chimeras in low-input samples
CCC_detail_chim_corr %>%
  filter(Dil %in% c(1.0e+05, 1.0e+04, 6.0e+03), Chim_prop > 0)

# Plot the data
svg("Output/Plots/F03A_CCC-chimeras-prop_2024-09-30.svg", 
    width = 3.85, height = 2.5)
set.seed(1)
CCC_detail_chim_corr %>%
  # Replace zeros by a small number for plotting
  mutate(Chim_prop = replace(Chim_prop, Chim_prop == 0, 0.00001)) %>%
  ggplot(., aes(x = Dil, y = Chim_prop)) +
  plot_theme +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept = 2e-5, linetype = "dashed", size = 0.4) +
  annotate("text", x = 3e7, y = 2.5e-5, hjust = 0, vjust = 0, label = "B.D.",
           size = 3) +
  geom_smooth(method = "lm", se = TRUE, colour = "#2D82C4", fill = "#69ABDD") + 
  geom_jitter(aes(shape = Broad_type), colour = "#174264", width = 0.12, 
              height = 0.12, size = 1.8) +
  xlab("Bacterial input cells [log10]") + ylab("Chimera rel. abund. [log10]") +
  annotate("text", x = 6000, y = 0.09, label = "r[S] == 0.88", 
           parse = TRUE, hjust = 0, vjust = 1, size = 3) +
  annotate("text", x = 6000, y = 10^-1.5, label = "p < 0.0001", 
           hjust = 0, vjust = 1, size = 3) +
  ggtitle("  Chimera formation by dilution") +
  theme(axis.title = element_text(size = rel(0.9)),
        plot.caption = element_text(hjust = 0), plot.title.position = "plot") +
  scale_shape_manual("Mock community", values = c(16, 15, 17))
dev.off()

################################################################################
#
# Unknowns
#
################################################################################

# ------------------------------------------------------------------------------
# Which samples contain which types of unknown taxa?
# ------------------------------------------------------------------------------

# Extract the main unknown ASVs in mock samples and negative controls
CCC_detail_unk <- CCC_detail %>%
  filter(Category == "Unknown") %>% 
  # Add mean relative abundance in skin samples
  mutate(Skin_abund = rowMeans(across(S.ID$skin))) %>%
  select(ASV_ID, Genus, Skin_abund, all_of(S.ID$notskin)) %>% 
  # Transform data to long format
  melt(., id.vars = c("ASV_ID", "Genus", "Skin_abund"),
       variable.name = "Sample_ID", value.name = "Rel_abund") %>% 
  # Add span for filtering
  filter(Rel_abund != 0) %>%
  group_by(ASV_ID) %>%
  mutate(span = n()) %>%
  ungroup() %>% 
  # Order samples and ASVs (for set.seed())
  arrange(Sample_ID, Rel_abund, ASV_ID) %>%
  as.data.frame()

# Classify mean relative abundance in skin samples for visualisation
CCC_detail_unk <- CCC_detail_unk %>%
  mutate(skin.cut = cut(Skin_abund, right = T,
                        breaks = c(-1, 0, 0.0001, 0.001, 0.01, 0.1, 1))) %>%
  mutate(skin.cut = recode_factor(skin.cut, `(-1,0]` = "0"))

# Create dataset in wide format for clustering
Cluster_unk <- CCC_detail_unk %>%
  # Only keep ASVs with span > 2
  filter(span > 2) %>%
  select(Sample_ID, ASV_ID, Rel_abund) %>% 
  # Order samples and ASVs (for set.seed())
  pivot_wider(names_from = ASV_ID, values_from = Rel_abund) %>%
  column_to_rownames("Sample_ID")
# Replace zeros with a low value to allow log transformation
format(min(unlist(Cluster_unk), na.rm = T), scientific = F)
Cluster_unk[is.na(Cluster_unk)] <- 0.00001

# Determine number of clusters with kmeans (k = 1 to 10)
set.seed(1)
tot_withinss <- map_dbl(1:10,  function(k){
  model <- kmeans(x = log10(t(Cluster_unk)), centers = k)
  model$tot.withinss})
# Plot the elbow plot
ggplot(data.frame(k = 1:10, tot_withinss = tot_withinss), 
       aes(x = k, y = tot_withinss)) +
  geom_line() + geom_point() + scale_x_continuous(breaks = 1:10)
# Plot the heatmap with 4 clusters
set.seed(1)
Heatmap(as.matrix(log10(Cluster_unk)), column_split = 4)
dev.off()

# Extract 4 clusters
set.seed(1)
ASV_clust_ord <- column_order(Heatmap(as.matrix(log10(Cluster_unk)), 
                                      column_split = 4))
# Connect order with ASV number and group name
Group_df <- data.frame(
  as.data.frame(melt(ASV_clust_ord, value.name = "ID")),
  ASV_ID = colnames(Cluster_unk)[melt(ASV_clust_ord, value.name = "ID")$ID]) %>%
  # Create group names
  mutate(L1 = case_when(L1 == 4 ~ "Skin cross-contamination",
                        L1 == 3 ~ "ZymoResearch buffer",
                        L1 == 2 ~ "Qiagen buffer",
                        L1 == 1 ~ "Inconclusive"))

# Merge relative abundance data with cluster information and sample metadata
CCC_detail_unk <- CCC_detail_unk %>%
  merge(., Group_df, by = "ASV_ID", all = TRUE) %>% 
  # Add NA for samples with zero rel. abundance of contaminating ASVs
  complete(Sample_ID, L1) %>% 
  merge(., Metafile[Metafile$Broad_type != "Skin", 
                    c("Sample_ID", "Prot", "Broad_type", "Dil")], 
        by = "Sample_ID", all = TRUE) %>%
  # Keep ASVs more or less sorted by relative abundance for plotting
  arrange(ASV_ID)

# Separate protocol information
CCC_detail_unk <- CCC_detail_unk %>%
  mutate(Kit = gsub("_._.", "", Prot),
         Lysis = str_extract(Prot, "S|T"),
         Buffer = gsub("._._", "", Prot),# ) %>% View#,
         # Italicize genera
         ASV_ID2 = paste0("italic('", Genus, "')~(", ASV_ID, ")"),
         ASV_ID2 = gsub("italic\\(\\'NA\\'\\)", "'NA'", ASV_ID2), 
         # Create longer sample info except for DNA samples and PCR controls
         Sample_info = case_when(grepl("DNA", Kit) ~ paste0(Kit, " (", Broad_type, ")"),
                                 grepl("PCR", Kit) ~ Kit,
                                 TRUE ~ paste0(Dil, " (", Broad_type, ")")))
# Shorten sample info to create x axis labels
CCC_detail_unk <- CCC_detail_unk %>%
  mutate(Sample_info = gsub("Even mock)", "ev.)", Sample_info),
         Sample_info = gsub("Staggered mock", "st.", Sample_info),
         Sample_info = gsub("Spike-in", "sp.", Sample_info),
         Sample_info = gsub("cells ", "", Sample_info),
         Sample_info = gsub("^6\\*", "", Sample_info),
         Sample_info = gsub("Pipel.*ntrol\\)", "N. control", Sample_info))
# Create facet labels for DNA samples and PCR controls
CCC_detail_unk[grepl("DNA", CCC_detail_unk$Kit), "Buffer"] <- "D"
CCC_detail_unk[grepl("DNA", CCC_detail_unk$Kit), "Lysis"] <- " "
CCC_detail_unk[grepl("DNA", CCC_detail_unk$Kit), "Kit"] <- "  "
CCC_detail_unk[grepl("PCR", CCC_detail_unk$Kit), "Buffer"] <- "A"
CCC_detail_unk[grepl("PCR", CCC_detail_unk$Kit), "Lysis"] <- " "
CCC_detail_unk[grepl("PCR", CCC_detail_unk$Kit), "Kit"] <- "  "

# Create a named colour vector for abundance in ComplexHeatmap
Col_abund <- c("#F5CBD8", "#EA96B1", "#DE6088", "#C90040", "#91002E", "#4B0018")
names(Col_abund) <- levels(CCC_detail_unk$skin.cut)

# Plot the data
svg("Output/Plots/S04_CCC-contaminant_2024-09-30.svg", width = 8.5, height = 3.7)
CCC_detail_unk %>%
  filter(span > 2 | is.na(span)) %>%
  # Fix order of facets, taxa, and samples 
  mutate(L1 = factor(L1, levels = c("Skin cross-contamination", "Qiagen buffer", 
                                    "ZymoResearch buffer", "Inconclusive")),
         ASV_ID2 = case_when(ASV_ID2 == "'NA'~(NA)" & grepl("cross", L1) ~ 
                               "italic('Cutibacterium')~(ASV_0003)",
                             TRUE ~ ASV_ID2),
         ASV_ID2 = factor(ASV_ID2, levels = rev(unique(ASV_ID2))),
         Sample_info = factor(Sample_info, levels = c(
           "10^8 (ev.)", "10^8 (st.)", "10^6 (ev.)", "10^6 (st.)", "10^5 (ev.)", 
           "10^5 (sp.)", "10^4 (ev.)", "10^3 (sp.)", "N. control", "DNA_a (ev.)", 
           "DNA_b (ev.)", "DNA_a (st.)", "DNA_b (st.)", "PCR_a", "PCR_b"))) %>%
  # Remove remaining NA values from previous "complete" command
  filter(ASV_ID2 != "'NA'~(NA)") %>%
  ggplot(., aes(x = Sample_info, y = ASV_ID2, size = Rel_abund, 
                colour = skin.cut)) +
  geom_point() +
  scale_y_discrete(breaks = CCC_detail_unk$ASV_ID2, 
                   labels = parse(text = CCC_detail_unk$ASV_ID2)) +
  ggh4x::facet_nested(L1 ~ Buffer + Kit + Lysis, scales = "free", space = "free") +
  scale_size("Relative abundance", range = c(0.3, 1.8)) +
  xlab("Dilution/replicate (ordered by decreasing input cells)") + 
  ylab("Genus (ASV)") + 
  scale_colour_manual("Mean rel. ab. in skin samples", values = Col_abund) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5,
                                   size = rel(0.45)),
        axis.text.y = element_text(size = rel(0.6)),
        axis.title = element_text(size = rel(0.8)),
        legend.title = element_text(size = rel(0.8)),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(size = rel(0.7)),
        strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, 
                                    size = rel(0.7)))
dev.off() # 8*8 + 8*Neg + 4*DNA + 2*PCR = 78 Samples

# ------------------------------------------------------------------------------
# Correlation of unknown proportion & the number of input cells per sample
# ------------------------------------------------------------------------------

# Select contaminants of all mock samples
CCC_detail_unk_prop <- CCC_detail_unk %>% 
  # Exclude negative controls
  filter(Sample_ID %in% S.ID$mock.all) %>%
  # Classify unknowns into skin cross contaminants and other contaminants
  mutate(Cont_group = case_when(
    L1 == "Skin cross-contamination" ~ "Cross-contaminants",
    TRUE ~ "Other contaminants"))

# Summarise relative abundances for cross- and other contaminants per sample
CCC_detail_unk_prop <- CCC_detail_unk_prop %>%
  group_by(Cont_group, Sample_ID, Kit, Lysis, Buffer, Dil, Broad_type) %>%
  summarise(val_sum = sum(Rel_abund, na.rm = TRUE))

# Replace zeros for plotting
CCC_detail_unk_prop$val_sum[CCC_detail_unk_prop$val_sum == 0] <- 0.00001
# Replace x axis labels for plotting
CCC_detail_unk_prop$Dil <- as.character(CCC_detail_unk_prop$Dil)
CCC_detail_unk_prop$Dil[CCC_detail_unk_prop$Dil == "D"] <- "1.3e7"
CCC_detail_unk_prop$Dil <- gsub(" cells", "", CCC_detail_unk_prop$Dil)
CCC_detail_unk_prop$Dil <- gsub("0\\^", "e", CCC_detail_unk_prop$Dil)
CCC_detail_unk_prop$Dil <- gsub("6\\*1", "6", CCC_detail_unk_prop$Dil)
CCC_detail_unk_prop$Dil <- as.numeric(CCC_detail_unk_prop$Dil)

# Spearman correlation test
cor.test(~ val_sum + Dil, method = "spearman",
         data = filter(CCC_detail_unk_prop, Cont_group == "Cross-contaminants"))
cor.test(~ val_sum + Dil, method = "spearman", 
         data = filter(CCC_detail_unk_prop, Cont_group == "Other contaminants"))

# Plot the data
svg("Output/Plots/F03C_CCC-cont-prop_2024-09-30.svg", width = 4.05, height = 2.5)
set.seed(1)
CCC_detail_unk_prop %>%
  ggplot(., aes(x = Dil, y = val_sum, colour = Cont_group, fill = Cont_group)) + 
  geom_jitter(aes(shape = Broad_type), width = 0.12, 
              height = 0.12, size = 1.8) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept = 2e-5, linetype = "dashed", size = 0.4) +
  annotate("text", x = 2e4, y = 2.5e-5, hjust = 1, vjust = 0, label = "B.D.",
           size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_colour_manual(name = NULL, values = c("#4B0018", "#D94574")) +
  scale_fill_manual(name = NULL, values = c("#7D2E47", "#EA99B3")) +
  plot_theme +
  ggtitle("      Contamination by dilution") +
  xlab("Bacterial input cells [log10]") + ylab("Contaminant rel. abund. [log10]") +
  theme(axis.title = element_text(size = rel(0.9)),
        plot.caption = element_text(hjust = 0), plot.title.position = "plot",
        legend.margin = margin(c(1,5,1,5)),
  ) +
  scale_shape_manual("Mock community", values = c(16, 15, 17)) +
  guides(shape = guide_legend(order = 1), 
         colour = guide_legend(order = 2), fill = guide_legend(order = 2)) +
  annotate("text", x = 10^8, y = 0.43, label = "r[S] < -0.70", 
           parse = TRUE, hjust = 1, vjust = 1, size = 3) +
  annotate("text", x = 10^8, y = 0.12, label = "p < 0.0001", 
           hjust = 1, vjust = 1, size = 3)
dev.off() # 8*8 + 4*DNA = 68 *2 = 136 Samples

# ------------------------------------------------------------------------------
# Contaminant composition plot (stacked bar chart)
# ------------------------------------------------------------------------------

# Classify all unknown ASVs and summarise proportions per sample
CCC_detail_unk_stack <- CCC_detail_unk %>%
  # Classify ASVs with span <=2 as "inconclusive"
  mutate(ASV_group = case_when(is.na(L1) ~ "Inconclusive", TRUE ~ L1)) %>%
  # Keep only negative controls
  filter(!Dil %in% c("A", "D"), grepl("Neg", Broad_type)) %>% 
  # Summarise contaminants per group (cross-cont./buffer/inconcl) per sample
  group_by(ASV_group, Sample_ID, Kit, Lysis, Buffer, Dil, Broad_type) %>%
  summarise(Rel_abund = sum(Rel_abund, na.rm = TRUE))

# Text: proportion of explained contaminants in pipeline negative controls
CCC_detail_unk_stack %>%
  group_by(Sample_ID, Kit, Lysis, Buffer, Dil, Broad_type) %>%
  # Scale to 100% to account for expected ASVs in negative controls 
  mutate(value_scaled = Rel_abund/sum(Rel_abund)) %>% 
  # Group contaminants into conclusive and inconclusive origin
  mutate(Explained = case_when(ASV_group == "Inconclusive" ~ "Inconcl",
                               TRUE ~ "Concl")) %>%
  # Calculate proportion of concl./inconcl. reads per control sample
  group_by(Explained, Sample_ID) %>%
  summarise(value_scaled = sum(value_scaled)) %>%
  # Calculate median proportion of reads
  group_by(Explained) %>%
  summarize(med = median(value_scaled)) # 47.3 % from cross-cont & buffers

# Plot the data
svg("Output/Plots/F03D_CCC-cont-stack_2023-08-08.svg", width = 4.335, height = 3.245) #width = 4.35, height = 3.27
CCC_detail_unk_stack %>%
  # Order fill categories (contaminant categories)
  mutate(ASV_group = factor(
    ASV_group, 
    levels = sort(unique(CCC_detail_unk$L1))[c(2, 4, 1, 3)])) %>%
  ggplot(aes(x = interaction(Kit, Lysis, Buffer), y = Rel_abund, fill = ASV_group)) +
  geom_bar(position = "fill", stat = "identity") +
  ggh4x::facet_nested(. ~ Buffer + Kit + Lysis, scales = "free", 
                      space = "free", switch = "x") +
  plot_theme + 
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = rel(0.9)),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(margin = margin(0, 2, 0, 0)),
        plot.title.position = "plot",
        panel.spacing = unit(0, "lines"),
        plot.tag.position = c(0.58, 0.18),
        plot.tag = element_text(size = rel(0.9), lineheight = 1.42, hjust = 0),
        plot.title = element_text(margin = margin(0, 0, 20, 0))) +
  scale_y_continuous(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) +
  xlab("Extraction protocol") + ylab("Proportion of contaminants") +
  scale_fill_manual("Contaminant origin", 
                    values = c("#FF0A68", "#FFAE00", "#FDD9E4", "#A97585")) +
  ggtitle("Contamination by protocol (pipeline neg. controls)") +
  labs(tag = "Lysis condition \nExtraction kit \nBuffer")
dev.off()

rm(ASV_clust_ord, CCC_detail, CCC_detail_chim, CCC_detail_chim_corr, 
   CCC_detail_chim_even, CCC_detail_chim_log, CCC_detail_chim_plot, 
   CCC_detail_err, CCC_detail_unk, CCC_detail_unk_prop, CCC_detail_unk_stack,
   Cluster_unk, Group_df, ASV_order, tax_order, tax_order2, tax1, tax2,
   tot_withinss)
