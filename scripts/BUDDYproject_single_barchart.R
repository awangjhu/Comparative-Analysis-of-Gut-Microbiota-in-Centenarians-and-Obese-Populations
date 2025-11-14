library(tidyverse)

#change file paths as needed
results_dir <- "~/qbb2025/BUDDYproject/results"
fig_dir <- "~/qbb2025/BUDDYproject/figures"

#open the data
cent <- read.csv(results_dir, "centenarian_S_7G6_report.csv", header=FALSE)
young <- read.csv(results_dir, "young_S_7G39_report.csv", header=FALSE)


#add designated column headers
colnames(cent) <- c("percent", "reads", "clade_reads", "rank_code", "ncbi_id", "name")
colnames(young) <- c("percent", "reads", "clade_reads", "rank_code", "ncbi_id", "name")


#filter centenarian cohort output data by class taxa
cent_class <- cent %>%
  filter(rank_code == "C") %>%
  mutate(sample = "Centenarian",
         percent = as.numeric(percent))


#filter young cohort output data by class taxa
young_class <- young %>%
  filter(rank_code == "C") %>%
  mutate(sample = "Young",
         percent = as.numeric(percent))

#combine the two output datas grouped by class
#normalize to percent composition for each cohort
#group percent by into other category
combined_class <- bind_rows(cent_class, young_class) %>%
  group_by(sample) %>%
  mutate(percent = percent / sum(percent) * 100)%>%
  mutate(name = ifelse(percent < 1, "Other", name)) %>%
  group_by(sample, name) %>%
  summarise(percent = sum(percent), .groups = "drop")%>%
  mutate(name = str_trim(name))

#barchart to visualize percent composition differences by taxanomic class
plot <- ggplot(combined_class, aes(x = sample, y = percent, fill = name)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_x_discrete(expand = expansion(mult = c(0.3, 0.3))) +
  labs(
    title = "Class-level Taxonomic Composition",
    x = "",
    y = "Relative Abundance (%)",
    fill = "Taxonomic Class"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text( vjust = 1, face = "bold"),
    plot.title = element_text(face = "bold", size = 16),  # centers + enlarges title
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),  # smaller legend text
    legend.key.size = unit(0.4, "cm")     # smaller color boxes
  )

#save image to filepath established earlier
ggsave(
  file.path(fig_dir, "class_composition_single.png"), plot,
  width = 8, height = 6, dpi = 300
)
