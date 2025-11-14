library(tidyverse)

# Define function to read and clean a single Kraken2 report file
read_report <- function(file) {
  df <- read.table(file, header = FALSE, fill = TRUE, quote = "", comment.char = "", sep = "")
  
  # If fewer than 6 columns, pad with NAs
  if (ncol(df) < 6) {
    missing_cols <- 6 - ncol(df)
    df[, (ncol(df) + 1):6] <- NA
  }
  # If more than 6 columns, collapse extras into the last column
  if (ncol(df) > 6) {
    df <- df[, 1:6]
  }
  
  colnames(df) <- c("percent", "reads", "clade_reads", "rank_code", "ncbi_id", "name")
  
  df <- df %>%
    mutate(
      percent = as.numeric(percent),
      name = str_trim(name)
    ) %>%
    filter(!is.na(percent))
  
  return(df)
}

# Define directories
results_directory <- "~/qbb2025/BUDDYproject/results"
avg_dir <- "~/qbb2025/BUDDYproject/avg_results"

# Find all centenarian report files
cent_files <- list.files(
  results_directory,
  pattern = "^centenarian_.*_report\\.(csv|txt)$",
  full.names = TRUE
)

print(cent_files)

# Combine and normalize each report
cent_combined <- map_dfr(cent_files, function(f) {
  read_report(f) %>%
    mutate(
      replicate_id = tools::file_path_sans_ext(basename(f))
    ) %>%
    group_by(replicate_id) %>%
    mutate(percent = percent / sum(percent, na.rm = TRUE) * 100) %>%  # normalize within sample
    ungroup()
})

# Average across replicates
cent_avg <- cent_combined %>%
  group_by(rank_code, ncbi_id, name) %>%
  summarise(mean_percent = mean(percent, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_percent))

# Save averaged results
write.table(
  cent_avg,
  file = file.path(avg_dir, "centenarian_average.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


young_files <- list.files(
  results_directory,
  pattern = "^young_.*_report\\.(csv|txt)$",
  full.names = TRUE
)

print(young_files)

young_combined <- map_dfr(young_files, function(f) {
  read_report(f) %>%
    mutate(replicate_id = tools::file_path_sans_ext(basename(f))) %>%
    group_by(replicate_id) %>% 
    mutate(percent = percent / sum(percent, na.rm = TRUE) * 100) %>%  # normalize within each replicate
    ungroup()
})

young_avg <- young_combined %>%
  group_by(rank_code, ncbi_id, name) %>%
  summarise(mean_percent = mean(percent, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_percent))

write.table(
  young_avg,
  file = file.path(avg_dir, "young_average.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
