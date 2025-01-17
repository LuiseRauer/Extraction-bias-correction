################################################################################
#
# Zymo-Qiagen comparison mock analysis
# Common functions and colour schemes
#
################################################################################

# Load required packages
library(tidyverse)

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

# Reverse complement of a single DNA sequence
rev_comp <- function(x) {
  # Separate letters 
  x <- strsplit(x, split = "")[[1]]
  # Reverse
  x <- rev(x)
  # Complement
  x <- chartr("ACGT", "TGCA", x)
  # Collapse letters
  x <- paste0(x, collapse = "")
  return(x)
}

# Geometric mean of a vector
geom_mean <- function(x) {exp(mean(log(x)))}

# Maximum relative abundance per row
rowMax <- function(x) {
  apply(x, 1, function(y) max(y, na.rm = TRUE))
}

# Span per row (values > 0)
rowSpan <- function(x) {
  apply(x, 1, function(y) sum(y > 0, na.rm = TRUE))
}

# ------------------------------------------------------------------------------
# Expected data
# ------------------------------------------------------------------------------

# Define expected mock taxa and relative abundances
Expected_mock <- data.frame(
  Species = c("Pseudomonas_aeruginosa", "Escherichia_coli",
              "Salmonella_enterica", "Lactobacillus_fermentum",
              "Enterococcus_faecalis", "Staphylococcus_aureus",
              "Listeria_monocytogenes", "Bacillus_subtilis"),
  Abund_even_16S = c(4.2, 10.1, # P. aeruginosa, E. coli,
                     10.4, 18.4, # S. enterica, L. fermentum,
                     9.9, 15.5, # E. faecalis, S. aureus,
                     14.1, 17.4), # L. monocytogenes, B. subtilis
  Abund_log_16S = c(2.8, 0.069, # P. aeruginosa, E. coli,
                    0.07, 0.012, # S. enterica, L. fermentum,
                    0.00067, 0.0001, # E. faecalis, S. aureus,
                    95.9, 1.2)) # L. monocytogenes, B. subtilis

# Define expected spike-in taxa and relative abundances (needed in script 8)
Expected_spike <- data.frame(
  Species = c("Truepera_radiovictrix", "Imtechella_halotolerans",
              "Allobacillus_halotolerans"),
  Abund_16S = c(84.2, 12.9, 2.9)) # Truepera, Imtech., Allobacillus

# ------------------------------------------------------------------------------
# Font setting of bacteria
# ------------------------------------------------------------------------------

species_labels <- c(
  "Others" = "Non-mock taxa (LV > 4)",
  "Unknown" = "Unknown",
  "Bacillus_subtilis" = expression(italic("Bacillus subtilis")),
  "Enterococcus_faecalis" = expression(italic("Enterococcus faecalis")),
  "Escherichia_coli" = expression(italic("Escherichia coli")),
  "Lactobacillus_fermentum" = expression(italic("Lactobacillus fermentum")),
  "Listeria_monocytogenes" = expression(italic("Listeria monocytogenes")), 
  "Pseudomonas_aeruginosa" = expression(italic("Pseudomonas aeruginosa")),
  "Salmonella_enterica" = expression(italic("Salmonella enterica")),
  "Staphylococcus_aureus" = expression(italic("Staphylococcus aureus")),
  "Allobacillus_halotolerans" = expression(italic("Allobacillus halotolerans")), 
  "Imtechella_halotolerans" = expression(italic("Imtechella halotolerans")),
  "Truepera_radiovictrix" = expression(italic("Truepera radiovictrix")))

# ------------------------------------------------------------------------------
# General plot themes and colours
# ------------------------------------------------------------------------------

# Define a common plot theme
plot_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "grey60", fill = NA),
  strip.background = element_rect(
    colour = "grey60", fill = "grey92", linewidth = 0.5, linetype = "solid"),
  legend.key = element_blank())

# Define colours for mock and spike-in species
Col_spec <- setNames( # used in sample composition!
  c("#FADD00", # yellow, Bacillus
    "#FF5C7F", # pink, Enterococcus
    "#23A959", # turquoise, Escherichia
    "#F5A700", # gold, Lactobacillus
    "#DB0042", # red, Listeria
    "#6BDB96", # light turquoise, Pseudomonas
    "#006631", # dark green, Salmonella
    "#FF99AF", # rose, Staphylococcus
    "#AD8002", # gold, Allobacillus
    "#697300", # olive, Imtechella
    "#8C3B32", # rust, Truepera
    "grey60"), # grey, Unknown
  nm = c(sort(Expected_mock$Species), sort(Expected_spike$Species), "Others"))

# Define colours for the classification of sequence errors & chimeras
Col_class <- c("#005710", "#4B8D15", "#8EB61F", "#D1DE19", "#FFF833", # greens
               "#FFC31F", "#EC7F32", "#D73030", "#A00000", # reds
               "#2D82C4", # chimera blue
               "grey60") # unknown grey

# Define colours for number of ASVs
Col_ASV <- c("#000000", "#2A0F44", "#531D87", "#972FB1", "#C93D83", "#E85151", 
             "#F1835B", "#F6B765", "#FAD975")


# Define colours for dilutions
Col_dil <- setNames(
  c("#1E0C46", "#4359B1", "#16BCF3", "#73E8E4", "#5CDBAE"),
  nm = c("10^8 cells", "10^6 cells", "10^5 cells", "10^4 cells", "6*10^3 cells")) 

# Define colours for protocols
Col_prot <- setNames(c("#2E58B8", "#5DAADE", # dark blue, light blue
                       "#009C00", "#7EDE68", # dark green, light green
                       "#ba002b", "#ffbdc7", # dark red, rose
                       "#E66C1E", "#ffda0a", # orange, yellow
                       "#2F9D92", "#EC796F", # dark cyan, coral - kit
                       "#964A8B", "#D99E32", # plum, gold - lysis condition
                       "#804141", "#E8BD9B"), # brown, sand - buffer
                     nm = c("Q_S_q", "Q_S_z", "Q_T_q", "Q_T_z", 
                            "Z_S_q", "Z_S_z", "Z_T_q", "Z_T_z",
                            "Q", "Z", "S", "T", "q", "z"))

# R objects after running this script:
# "Col_ASV", "Col_class", "Col_dil", "Col_prot", "Col_spec", "Expected_mock", 
#   "Expected_spike", "geom_mean", "plot_theme", "rev_comp", "rowMax", 
#   "species_labels", "rowSpan"
