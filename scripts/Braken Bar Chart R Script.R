#!/usr/bin/env Rscript
# Class-level microbiome analysis for centenarian study
# Analyzes Bracken output comparing young vs elderly vs centenarian. Since n = 1, we are still working on
# importing the entire dataset from ENA

setwd("C:/Users/alanw/OneDrive - Johns Hopkins/Quant Bio Project")

library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)

# Load Data
# Load class-level relative abundance (individual samples)
class_rel <- read.table("results/merged/class_relative_abundance.tsv",
                        header = TRUE,
                        sep = "\t",
                        row.names = 1)

# Loading absolute abundance 
class_abs <- read.table("results/merged/class_abundance.tsv",
                        header = TRUE,
                        sep = "\t",
                        row.names = 1)

print("=== DATA LOADED ===")
print(paste("Total bacterial classes detected:", nrow(class_rel)))
print(paste("Samples:", paste(colnames(class_rel), collapse=", ")))

# Assigning sample groups from our metadata
# Sample metadata
sample_groups <- data.frame(
  sample = c("S_7G6", "S_7G38", "S_7G39"),
  age_group = c("Centenarian", "Elderly", "Young"),
  age = c(100, 75, 30)
)

print("\n=== SAMPLE GROUPS ===")
print(sample_groups)

#Comparison Table
comparison <- data.frame(
  taxon = rownames(class_abs),
  young = class_abs$S_7G39,
  elderly = class_abs$S_7G38,
  centenarian = class_abs$S_7G6
)

# Calculate fold changes
comparison$fold_change_cent_vs_young <- comparison$centenarian / comparison$young
comparison$fold_change_elderly_vs_young <- comparison$elderly / comparison$young

# Filter out low-abundance classes
min_abundance <- 100  # minimum reads in any group
comparison_filtered <- comparison %>%
  filter(young > min_abundance | elderly > min_abundance | centenarian > min_abundance)

print("\n=== ABUNDANCE FILTERING ===")
print(paste("Classes after filtering (>", min_abundance, "reads):", 
            nrow(comparison_filtered)))

# Identifying bacterial class
# Classes that decline from young to centenarian (>50% reduction)
declining <- comparison_filtered %>%
  filter(fold_change_cent_vs_young < 0.5) %>%
  arrange(fold_change_cent_vs_young)

print("\n=== DECLINING BACTERIAL CLASSES (>50% reduction in centenarians) ===")
if(nrow(declining) > 0) {
  print(declining[, c("taxon", "young", "elderly", "centenarian", "fold_change_cent_vs_young")])
} else {
  print("No classes showed >50% decline")
}

# Classes that increase
increasing <- comparison_filtered %>%
  filter(fold_change_cent_vs_young > 2.0) %>%
  arrange(desc(fold_change_cent_vs_young))

print("\n=== INCREASING BACTERIAL CLASSES (>2x in centenarians) ===")
if(nrow(increasing) > 0) {
  print(increasing[, c("taxon", "young", "elderly", "centenarian", "fold_change_cent_vs_young")])
} else {
  print("No classes showed >2-fold increase")
}


# Stacked Bar Chart
# Select top 15 most abundant classes across all samples
top_classes <- rownames(class_rel)[order(rowSums(class_rel), decreasing = TRUE)[1:15]]

# Data
class_rel_subset <- class_rel[top_classes, ]
class_rel_subset$taxon <- rownames(class_rel_subset)

# ggplot
class_long <- class_rel_subset %>%
  pivot_longer(cols = -taxon, names_to = "sample", values_to = "rel_abundance") %>%
  left_join(sample_groups, by = "sample") %>%
  mutate(age_group = factor(age_group, levels = c("Young", "Elderly", "Centenarian")))

# Color palette
n_taxa <- length(top_classes)
colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_taxa)

# Stacked Bar Chart
p1 <- ggplot(class_long, aes(x = age_group, y = rel_abundance, fill = taxon)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  scale_fill_manual(values = colors) +
  labs(
    title = "Class-Level Taxonomic Composition of Gut Microbiota",
    subtitle = "Sardinia Centenarian Study (n=1 per group)",
    x = "Age Group",
    y = "Relative Abundance (%)",
    fill = "Bacterial Class"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold")
  )

ggsave("results/Figure1_class_stacked_bar.pdf", p1, width = 10, height = 7)
ggsave("results/Figure1_class_stacked_bar.png", p1, width = 10, height = 7, dpi = 300)

print("\n=== FIGURE 1 SAVED ===")
print("results/Figure1_class_stacked_bar.pdf")
print("results/Figure1_class_stacked_bar.png")