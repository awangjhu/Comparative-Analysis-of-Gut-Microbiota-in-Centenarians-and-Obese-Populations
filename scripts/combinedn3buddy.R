library(tidyverse)

avg_dir <- "~/qbb2025/BUDDYproject/avg_results"

cent_avg <- read.table(file.path(avg_dir, "centenarian_average.txt"), header = TRUE, sep = "\t")
young_avg <- read.table(file.path(avg_dir, "young_average.txt"), header = TRUE, sep = "\t")


cent_avg <- cent_avg %>%
  mutate(group = "Centenarian")

young_avg <- young_avg %>%
  mutate(group = "Young")


combined_avg <- bind_rows(cent_avg, young_avg)


combined_class <- combined_avg %>%
  filter(rank_code == "C") %>%
  group_by(group) %>%
  mutate(mean_percent = mean_percent / sum(mean_percent) * 100) %>%
  mutate(name = ifelse(mean_percent < 1, "Other", name)) %>%
  group_by(group, name) %>%
  summarise(mean_percent = sum(mean_percent), .groups = "drop") %>%
  mutate(name = str_trim(name))


ggplot(combined_class, aes(x = group, y = mean_percent, fill = name)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_x_discrete(expand = expansion(mult = c(0.3, 0.3))) +
  labs(
    title = "Class-level Taxonomic Composition 
                        (Averaged)",
    x = "",
    y = "Relative Abundance (%)",
    fill = "Taxonomic Class"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(vjust = 1, face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

