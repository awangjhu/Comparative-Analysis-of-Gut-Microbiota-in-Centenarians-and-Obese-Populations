library(tidyverse)
library(phyloseq)
library(ANCOMBC)


results_dir <- "~/qbb2025/BUDDYproject/brckn_results"
fig_dir <- "~/qbb2025/BUDDYproject/figures"

bracken_files <- list.files(
  bracken_dir,
  pattern = "bracken\\.report\\.txt$",
  full.names = TRUE
)


read_bracken <- function(file) {
  
  df <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  
  # Bracken standard format:
  # percent | reads | extra | rank | taxid | name
  colnames(df) <- c("percent", "reads", "extra", "rank", "taxid", "name")
  
  df %>%
    filter(rank == "C") %>%    # Class level only
    mutate(
      reads = as.numeric(reads),
      name = str_trim(name),
      sample = tools::file_path_sans_ext(basename(file))
    ) %>%
    select(sample, name, reads)
}
class_counts <- map_dfr(bracken_files, read_bracken)

class_matrix <- class_counts %>%
  pivot_wider(
    names_from = sample,
    values_from = reads,   # Bracken corrected counts
    values_fill = 0
  )



samples <- colnames(class_matrix)[-1]  # exclude "name"



metadata <- tibble(
  sample = samples,
  Group = ifelse(grepl("G6", sample, ignore.case = TRUE),  # example pattern
                 "Centenarian",
                 "Young")
)

metadata


otu <- as.data.frame(class_matrix)
rownames(otu) <- otu$name
otu <- otu[, -1]
otu <- as.matrix(otu)

OTU <- otu_table(otu, taxa_are_rows = TRUE)

TAX <- tax_table(matrix(
  rownames(otu),
  ncol = 1,
  dimnames = list(rownames(otu), "Class")
))

META <- sample_data(metadata)
rownames(META) <- META$sample

phy <- phyloseq(OTU, TAX, META)

phy


ancom_out <- ancombc2(
  data = phy,
  assay_name = "counts",
  fix_formula = "Group",
  rand_formula = NULL,
  group = "Group",       
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  struc_zero = TRUE
)

res <- ancom_out$res
head(res)

sig <- res %>%
  filter(q_GroupYoung < 0.05)

sig

volcano_data <- res %>%
  select(
    Class = taxon,
    logFC = lfc_GroupYoung,
    pval = p_GroupYoung,
    padj = q_GroupYoung
  )



plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(padj))) +
  geom_point(aes(color = logFC > 0), size = 3) +
  scale_y_continuous(limits = c(0, max(-log10(volcano_data$padj)) + 0.2)) +
  theme_minimal() +
  labs(
    title = "BrackenANCOM-BC2 Differential Abundance (Class Level)",
    x = "Log Fold Change (Young vs Centenarian)",
    y = "-log10(FDR)"
  ) +
  scale_color_manual(values = c("gray", "red"))

ggsave(
  file.path(fig_dir, "brckn_ancom_volcano.png"), plot,
  width = 8, height = 6, dpi = 300
)
