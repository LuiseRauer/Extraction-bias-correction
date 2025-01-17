################################################################################
#
# Bias analysis in skin samples
#
################################################################################

# Load required libraries
library(tidyverse)
library(reshape2)
library(openxlsx)
library(vegan)
library(ggh4x)
library(cowplot)

# Define species labels for skin taxa
species_labels2 <- c(
  "Others" = "Others",
  "Aerococcus" = expression(italic("Aerococcus")),
  "Blautia" = expression(italic("Blautia")),
  "Lawsonella" = expression(italic("Lawsonella")),
  "Staphylococcus" = expression(italic("Staphylococcus")),
  "Streptococcus" = expression(italic("Streptococcus")),
  "Corynebacterium" = expression(italic("Corynebacterium")), 
  "Cutibacterium" = expression(italic("Cutibacterium")),
  "Faecalibacterium" = expression(italic("Faecalibacterium")),
  "Acinetobacter" = expression(italic("Acinetobacter")), 
  "Paracoccus" = expression(italic("Paracoccus")),
  "Lactobacillus" = expression(italic("Lactobacillus")),
  "Sphingomonas" = expression(italic("Sphingomonas")))

# Define colours for skin genera
Col_spec2 <- Col_spec[grepl("Others|Staphy", names(Col_spec))]
names(Col_spec2) <- gsub("_.*", "", names(Col_spec2))

Col_spec2 <- c(Col_spec2, setNames(
  c("#FF0040", "#CD5555", "pink", "red3",
    "#F5DC00", "#F5B400", "gold3", "#F5A700",
    "#5CACEE", "#53DF82", "forestgreen"),
  nm = c("Lawsonella", "Streptococcus", "Aerococcus", "Blautia", 
         "Cutibacterium", "Faecalibacterium", "Corynebacterium", "Lactobacillus",
         "Paracoccus", "Acinetobacter", "Sphingomonas"))) # Staph already included

# ------------------------------------------------------------------------------
# Stacked bar chart before correction
# ------------------------------------------------------------------------------

# Check number of ASVs per skin sample
apply(seqtab.final[, S.ID$skin[c(1:6, 8:14, 16)]], 2, function(x) sum(x > 0))
# 5 - 11 - 106 - ... - 390 - 444 - 650

S.ID_skinsub <- Metafile %>%
  filter(Sample_ID %in% S.ID$skin, !grepl("Z_T", Prot)) %>% 
  select(Sample_ID) %>% unlist() %>% unname()

# Check annotation of buffer contaminants
seqtab.final.rel[grepl("0050|0034|0033|0042|0047|0056|0032|0043", 
                       seqtab.final.rel$ASV_ID), c(S.ID$skin, "LV_10", "Genus")]

# ASV table without contaminants, on required taxon level
seqtab.rel_skin <- seqtab.final.rel %>% 
  # Remove buffer contaminants and spike-in taxa
  filter(!grepl("0050|0034|0033|0042|0047|0056|0032|0043", ASV_ID)) %>%
  filter(!grepl("Methylobacterium|Brucella|Alcaligenes|Paraburkholder", Genus)) %>% 
  filter(!grepl("Nitrospirillum|Herbaspirillum|Aquabacterium", Genus)) %>% 
  filter(!grepl("Truepera|Imtechella|Allobacillus", Genus)) %>% # only 3263 remain!
  # Summarise counts by genus
  group_by(Genus) %>%
  summarise_at(.vars = S.ID_skinsub, .funs = sum) %>%
  as.data.frame()

# Top 10 skin genera
top_10_skin <- seqtab.rel_skin %>%
  mutate(meanVar = rowMeans(across(all_of(S.ID_skinsub)), na.rm = TRUE),
         span = rowSpan(across(all_of(S.ID_skinsub)))) %>%
  #filter(span >= 8) %>%
  arrange(desc(meanVar)) %>% #View
  top_n(11, meanVar) %>% # change between 22, 12
  select(Genus) %>% unlist(use.names = FALSE)
# Remove NA
top_10_skin <- na.omit(top_10_skin)

################################################################################
#
# Bias correction
#
################################################################################

# Read gram and shape information for skin taxa 
Skin_stain <- read.xlsx(paste0(path_proj, "/Extraction-bias-correction/Input/Skin-taxa_2024-09-04.xlsx"))

Bias_even_106 # Factor from 10^6 even mock sample

# Calculate correction of skin samples
data_skin_corr <- seqtab.final.rel[, c(S.ID_skinsub, "Genus")] %>% 
  mutate(Others = case_when(Genus %in% top_10_skin ~ Genus,
                           TRUE ~ "Others")) %>%
  group_by(Others) %>% 
  summarise_at(.vars = S.ID_skinsub, .funs = sum) %>% 
  melt(.) %>% 
  # Merge with skin data info
  merge(., Skin_stain, by.x = "Others", by.y = "Genus", all.x = TRUE) %>% 
  # Merge with protocol information?
  merge(., Metafile[, c("Sample_ID", "Prot", "Dil")], by.x = "variable", 
        by.y = "Sample_ID", all.x = TRUE) %>% #View
  merge(., Bias_even_106, by.x = c("Group", "Prot"),
        by.y = c("Group", "Prot"), all.x = TRUE) %>%
  #filter(Others != "Others", Group != "coccoid, gram-neg.") %>%
  mutate(value_corr = case_when(is.na(estimate_gm) ~ value, 
                                TRUE ~ value / estimate_gm)) %>% 
  # scale corrected value to 100%
  group_by(variable) %>%
  mutate(value_corr_scaled = value_corr/sum(value_corr),
         value_scaled = value/sum(value))

data_skin_corr1 <- data_skin_corr %>%
  select(Group, Prot, variable, Others, Dil, value_scaled, value_corr_scaled) %>%
  melt(id.vars = c("Group", "Prot", "variable", "Others", "Dil"), 
       variable.name = "Type") %>% 
  mutate(Type = case_when(Type == "value_scaled" ~ "Original", TRUE ~ "Corrected"),
         Type = factor(Type, levels = c("Original", "Corrected"))) %>%
  mutate(Others = factor(Others, levels = c("Others", sort(unique(.$Others))[sort(unique(.$Others)) != "Others"]))) %>%
  mutate(Dil_prot = paste0(Dil, ".", Prot), 
         Kit = gsub("_._.", "", Prot),
         Lysis = str_extract(Prot, "S|T"),
         Buffer = gsub("._._", "", Prot))

# Subject 1
pp1 <- data_skin_corr1 %>%
  filter(Dil == "Subject 1") %>%
  ggplot(., aes(x = Type, y = value, fill = Others)) +
  # Fill or stacked bar chart
  geom_bar(stat = "identity", position = "fill") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.35),
        legend.position = "none",
        axis.title.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.x = element_text(colour = "transparent")) +
  ylab("Relative abundance") + xlab("Status") +
  facet_nested(. ~ Dil + Kit + Lysis + Buffer, scales = "free_x") +
  scale_y_continuous(expand = c(0, 0), sec.axis = dup_axis()) + 
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual("Genus", values = Col_spec2) 

# Subject 2
pp2 <- data_skin_corr1 %>%
  filter(Dil == "Subject 2") %>%
  ggplot(., aes(x = Type, y = value, fill = Others)) +
  # Fill or stacked bar chart
  geom_bar(stat = "identity", position = "fill") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.35),
        legend.position = "none",
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_text(colour = "transparent"),
        axis.title.y.right = element_blank(),
        axis.text.y.right = element_text(colour = "white"),
        axis.ticks.y.right = element_line(colour = "transparent"),
        axis.title.x = element_text(colour = "transparent")) +
  ylab("Relative abundance") + xlab("Status") +
  facet_nested(. ~ Dil + Kit + Lysis + Buffer, scales = "free_x") +
  scale_y_continuous(expand = c(0, 0), sec.axis = dup_axis()) + 
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual("Genus", values = Col_spec2) 

# Legend
for (i in c("rods, gram neg.", "rods, gram pos.", "coccoid, gram pos.", 
            "coccoid, gram neg.", "Other")) {
  p <- data_skin_corr1 %>%
    mutate(Group = case_when(is.na(Group) ~ "Other", TRUE ~ Group)) %>%
    filter(Group == i) %>%
    ggplot(., aes(x = Type, y = value, fill = Others)) +
    # Fill or stacked bar chart
    geom_bar(stat = "identity", position = "fill") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.35),
          legend.text.align = 0,
          legend.title = element_text(size = rel(0.9)),
          legend.margin = margin(0.15, 0, 0.5, 0, unit = "cm"),
          legend.justification = c(0, 1))
  
  if(i == "rods, gram neg.") {p1 <- p + 
    scale_fill_manual("Rods, gram negative", values = Col_spec2, 
                      labels = species_labels2)}
  if(i == "rods, gram pos.") {p2 <- p +
    scale_fill_manual("Rods, gram positive", values = Col_spec2, 
                      labels = species_labels2)}
  if(i == "coccoid, gram pos.") {p3 <- p +
    scale_fill_manual("Ovoid, gram positive", values = Col_spec2, 
                      labels = species_labels2)}
  if(i == "coccoid, gram neg.") {p4 <- p +
    scale_fill_manual("Ovoid, gram negative", values = Col_spec2, 
                      labels = species_labels2)}
  if(i == "Other") {p5 <- p +
    scale_fill_manual("Not determined", values = Col_spec2, 
                      labels = species_labels2)}
}
# Build the legend
p_leg <- plot_grid(as_ggplot(get_legend(p3)), as_ggplot(get_legend(p2)),
                   as_ggplot(get_legend(p1)), as_ggplot(get_legend(p4)),
                   as_ggplot(get_legend(p5)),
                   ncol = 1, align = "hv", #rel_heights = c(1.33, 0.92, 0.5, 0.5, 0.5)
                   rel_heights = c(1.19, 1, 0.72, 0.5, 0.5, 0.5) #c(1.18, 0.73, 0.5, 0.5, 0.5)
                   )

# Build the plot: 
svg("Output/Plots/S08A_Skin-barchart.svg", width = 8, height = 4.2)
plot_grid(pp1, pp2, plot_grid(p_leg, NULL, nrow = 2, rel_heights = c(14, 1)), nrow = 1) +
  labs(tag = "Status") +
  theme(plot.tag.position = c(0.343, 0.05),
        plot.tag = element_text(size = 11, face = "plain", lineheight = 1.42, 
                                hjust = 0.5))
dev.off()

rm(pp1, pp2, p, p1, p2, p3, p4, p5, i)

################################################################################
#
# Change in relative abundance in skin samples - top three species
#
################################################################################

# ------------------------------------------------------------------------------
# Stacked bar chart before correction
# ------------------------------------------------------------------------------

# Top 3 skin species
top_3_skin <- seqtab.final.rel %>% 
  # Remove buffer contaminants and spike-in taxa
  filter(!grepl("0050|0034|0033|0042|0047|0056|0032|0043", ASV_ID)) %>%
  filter(!grepl("Methylobacterium|Brucella|Alcaligenes|Paraburkholder", Genus)) %>% 
  filter(!grepl("Nitrospirillum|Herbaspirillum|Aquabacterium", Genus)) %>% 
  filter(!grepl("Truepera|Imtechella|Allobacillus", Genus)) %>% # only 3263 remain!
  # Summarise counts by genus
  group_by(Species) %>%
  summarise_at(.vars = S.ID$skin[c(1:6, 8:14, 16)], .funs = sum) %>%
  as.data.frame() %>%
  mutate(meanVar = rowMeans(.[S.ID$skin[c(1:6, 8:14, 16)]], na.rm = TRUE),
         span = rowSpan(across(S.ID$skin[c(1:6, 8:14, 16)]))) %>%
  #filter(span >= 8) %>%
  arrange(desc(meanVar)) %>% 
  top_n(4, meanVar) %>% # change between 22, 12
  select(Species) %>% unlist(use.names = FALSE)
# Remove NA
top_3_skin <- na.omit(top_3_skin)
sort(top_3_skin)

################################################################################
#
# Bias correction
#
################################################################################

Bias_even_106 # Factor from 10^6 even mock sample

# Calculate correction of skin samples
data_skin_corr <- seqtab.final.rel[, c(S.ID_skinsub, "Species")] %>% 
  mutate(Others = case_when(Species %in% top_3_skin ~ Species,
                            TRUE ~ "Others"),
         Others_spec = gsub("_.*$", "", Others)) %>%
  group_by(Others, Others_spec) %>% 
  summarise_at(.vars = S.ID_skinsub, .funs = sum) %>% 
  melt(.) %>% 
  # Merge with skin data info
  merge(., Skin_stain, by.x = "Others_spec", by.y = "Genus", all.x = TRUE) %>% 
  # Merge with protocol information?
  merge(., Metafile[, c("Sample_ID", "Prot", "Dil")], by.x = "variable", 
        by.y = "Sample_ID", all.x = TRUE) %>% #View
  merge(., Bias_even_106, by.x = c("Group", "Prot"),
        by.y = c("Group", "Prot"), all.x = TRUE) %>%
  #filter(Others != "Others", Group != "coccoid, gram-neg.") %>%
  mutate(value_corr = case_when(is.na(estimate_gm) ~ value, 
                                TRUE ~ value / estimate_gm)) %>% 
  # scale corrected value to 100%
  group_by(variable) %>%
  mutate(value_corr_scaled = value_corr/sum(value_corr),
         value_scaled = value/sum(value))

################################################################################
#
# Corrected and not corrected data in one plot - rel. abundance plot with arrows
#
################################################################################

# Subject 1
p1 <- data_skin_corr %>%
  filter(Others != "Others",
         Dil == "Subject 1") %>%
  select(Group, Prot, variable, Others, Dil, value_scaled, value_corr_scaled) %>%
  mutate(Others = sub("^(([^_]*_){1}[^_]*).*", "\\1", Others),
         Others = gsub("_", "\n", Others),
         Others = factor(Others, levels = c("Cutibacterium\nacnes", 
                                            "Corynebacterium\ntuberculostearicum",
                                            "Staphylococcus\nepidermidis")),
         Group = factor(Group, levels = c("rods, gram pos.", "coccoid, gram pos.")),
         Kit = gsub("_._.", "", Prot),
         Lysis = str_extract(Prot, "S|T"),
         Buffer = gsub("._._", "", Prot)
  ) %>% 
  mutate(x = "1", y = "2") %>%
  filter(value_scaled != 0) %>% 
  ggplot(aes(x = x, xend = y, y = value_scaled, yend = value_corr_scaled, 
             colour = Others)) +
  geom_segment(arrow = arrow(length = unit(0.14, "cm")), lineend = "round", 
               linejoin = "round", size = 1.2) +
  geom_point(size = 1.5) +
  plot_theme +
  ylab("Relative abundance") + xlab("") +
  scale_x_discrete(labels = c("Original", "Corrected")) +
  scale_colour_manual(values = unname(Col_spec2[c(7, 9, 1)])) +
  facet_nested(Group + Others ~ Dil + Kit + Lysis + Buffer, scales = "free", 
               space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none") +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(0, 0.6)),
    scale_y_continuous(limits = c(0, 0.065)),
    scale_y_continuous(limits = c(0, 0.065))))

# Subject 2
p2 <- data_skin_corr %>%
  mutate(Group = case_when(Group == "coccoid, gram pos." ~ "ovoid, gram pos.",
                           TRUE ~ Group)) %>% 
  filter(Others != "Others",
         Dil == "Subject 2") %>%
  select(Group, Prot, variable, Others, Dil, value_scaled, value_corr_scaled) %>%
  mutate(Others = sub("^(([^_]*_){1}[^_]*).*", "\\1", Others),
         Others = gsub("_", "\n", Others),
         Others = factor(Others, levels = c("Cutibacterium\nacnes", 
                                            "Corynebacterium\ntuberculostearicum",
                                            "Staphylococcus\nepidermidis")),
         Group = factor(Group, levels = c("rods, gram pos.", "ovoid, gram pos.")),
         Kit = gsub("_._.", "", Prot),
         Lysis = str_extract(Prot, "S|T"),
         Buffer = gsub("._._", "", Prot)) %>% 
  mutate(x = "1", y = "2") %>%
  filter(value_scaled != 0) %>% 
  ggplot(aes(x = x, xend = y, y = value_scaled, yend = value_corr_scaled, 
             colour = Others)) +
  geom_segment(arrow = arrow(length = unit(0.14, "cm")), lineend = "round", 
               linejoin = "round", size = 1.2) +
  geom_point(size = 1.5) +
  plot_theme +
  ylab("Relative abundance") + xlab("") +
  scale_x_discrete(labels = c("Original", "Corrected")) +
  scale_colour_manual(values = unname(Col_spec2[c(7, 9, 1)])) +
  facet_nested(Group + Others ~ Dil + Kit + Lysis + Buffer, scales = "free", space = "free_x",
               strip = strip_nested(
                 text_y = elem_list_text(face = c(
                   "plain", # rods
                   "plain", # cocci
                   "italic", 
                   "italic",
                   "italic")))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        strip.text.y = element_text(angle = 0),
        axis.title.y = element_text(colour = "transparent", size = rel(1.25)),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none") +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(0, 0.6)),
    scale_y_continuous(limits = c(0, 0.065)),
    scale_y_continuous(limits = c(0, 0.065))))

# Join the plots
svg("Output/Plots/S08B_Skin-change.svg", width = 7.39, height = 3.7)
plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1.81)) +
  labs(tag = "          Status") +
  theme(plot.tag.position = c(0.343, 0.05),
        plot.tag = element_text(size = 11, face = "plain", lineheight = 1.42, 
                                hjust = 0.5))
dev.off()

# Table of correction factors
Bias_even_106 %>%
  mutate(estimate_gm = round(1/estimate_gm, 2)) %>%
  pivot_wider(names_from = "Prot", values_from = "estimate_gm") %>% 
  as.data.frame() %>% t() %>% 
  View

