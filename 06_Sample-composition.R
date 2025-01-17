################################################################################
#
# Zymo-Qiagen comparison mock analysis
# Sample composition plots
#
################################################################################

# Load required packages
library(tidyverse)
library(reshape2)

################################################################################
#
# Taxonomy plots
#
################################################################################

# Define a flexible function to create taxonomy plots
tax_plot_function <- function(
  data, # a seqtab
  tax_col, # column used for taxonomic annotation
  plot_op, # geom_bar position: "fill" or "stack"
  grid, # variable for facet_grid
  meta, # metadata
  colours # fill colours
  ) { 
  data %>%
    # Group by taxonomy variable & summarise counts
    group_by_at(tax_col) %>%
    summarise_if(is.numeric, sum) %>%
    # Create long format
    melt(., id.vars = tax_col) %>% 
    # Merge with metadata
    merge(., meta, by.x = "variable", by.y = "Sample_ID", all.x = TRUE) %>% 
    rename(grid = all_of(grid)) %>%
    mutate(variable = factor(variable, 
                             levels = arrange(Metafile, Prot)$Sample_ID)) %>%
    # Plotting
    ggplot(., aes_string(x = "Prot", y = "value", fill = tax_col)) +
    # Fill or stacked bar chart
    geom_bar(stat = "identity", position = plot_op) +
    plot_theme +
    # Facet
    ggh4x::facet_nested(. ~ Broad_type + eval(grid), 
                        scales = "free_x", space = "free") +
    # Rotate x axis labels
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                     size = 5),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_y_continuous(expand = c(0, 0)) + 
    ylab("Relative abundance") +
    guides(fill = guide_legend(ncol = 1))
}

# ------------------------------------------------------------------------------
# Taxonomy plot for sequence classification (sequence errors, chimeras, unknown)
# ------------------------------------------------------------------------------

# More precise classification of each sequence per sample
CCC_tax <- CCC_res %>% 
  mutate(Cat_LV = factor(case_when(
    Category == "Exact_match" ~ "LV = 0 (exact match)",
    Category == "Amplification_error" ~ paste0("LV = ", Min_LV_dist),
    Category == "Mock_chimera" ~ "LV > 8, mock chimera",
    Category == "Unknown" ~ "LV > 8, unclassified"), 
    levels = c("LV = 0 (exact match)", paste0("LV = ", 1:8), 
               "LV > 8, mock chimera", "LV > 8, unclassified"))) %>%
  group_by(Cat_LV) %>%
  summarise_at(.vars = sample.names, .funs = sum)

svg("Output/Plots/F02A_CCC-taxonomy_2023-06-29.svg", width = 10.48, height = 4)
tax_plot_function(
  CCC_tax, "Cat_LV", "fill", "Dil", Metafile %>% 
    mutate(Dil = factor(gsub(" cells", "", Dil), levels = 
                          gsub(" cells", "", levels(Metafile$Dil))))) +
  scale_fill_manual("Distance to exp. mock sequence", values = Col_class) 
dev.off()

# ------------------------------------------------------------------------------
# Text information
# ------------------------------------------------------------------------------

# Median proportion of exact matches in staggered mock with 10^6 cells
CCC_tax %>% 
  # Calculate relative abundances
  mutate(across(S.ID$all, ~ .x/sum(.x))) %>%
  filter(Cat_LV == "LV = 0 (exact match)") %>%
  # Select all staggered mock 10^6 samples
  select(S.ID$m.l.dil2) %>% 
  t(.) %>% median(.) # 99.7%

# Range of proportions of sequence errors in even mock over all dilutions
CCC_tax %>% 
  # Calculate relative abundances
  mutate(across(S.ID$all, ~ .x/sum(.x))) %>%
  filter(Cat_LV %in% paste0("LV = ", 1:8)) %>%
  # Select all even cell and DNA mock samples
  select(union(S.ID$mock.even, S.ID$DNA.even)) %>%
  colSums() %>% sort(.) # Range 6.2% - 31.6%

# Median proportion of sequence errors in even mock over all dilutions
CCC_tax %>% 
  # Calculate relative abundances
  mutate(across(S.ID$all, ~ .x/sum(.x))) %>%
  filter(Cat_LV %in% paste0("LV = ", 1:8)) %>%
  # Select all even cell and DNA mock samples
  select(union(S.ID$mock.even, S.ID$DNA.even)) %>%
  colSums() %>% sort(.) %>%  median(.) # Median 22.3%

# ------------------------------------------------------------------------------
# Taxonomy plot with final taxonomy with LV_4
# ------------------------------------------------------------------------------

# Replace NA by Others in LV_4 for plotting
seqtab.final <- seqtab.final %>%
  mutate(LV_4_Others = case_when(!is.na(LV_4) ~ LV_4, TRUE ~ "Others"),
         LV_4_Others = factor(LV_4_Others, levels = c(
           unique(Table_16S_279$ASV.ID), "Others")))

# LV_4, fill:
svg("Output/Plots/F02B_Taxonomy_LV4_2024-09-30.svg", width = 10, height = 4)
tax_plot_function(
  seqtab.final, "LV_4_Others", "fill", "Dil", Metafile %>% 
    mutate(Dil = factor(gsub(" cells", "", Dil), levels = 
                          gsub(" cells", "", levels(Metafile$Dil)))), Col_spec) +
  scale_fill_manual("Taxonomy", values = Col_spec, labels = species_labels) +
  theme(legend.text.align = 0)
dev.off()

rm(CCC_tax, tax_plot_function)
