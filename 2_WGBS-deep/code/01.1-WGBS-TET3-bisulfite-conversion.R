library(tidyverse)

# Chloroplast -----

# Load the chloroplast data
chloroplast_data <- read_tsv(
  "/path/to/project/WGBS-selected/data/4_methylation/chloroplast_data.tsv",
  show_col_types = FALSE,
  col_names = c("chromosome", "start", "end", "context", "methyation_level", "methylated_count", 
                "unmethylated_count", "sample")
)

# Calculate conversion rate for each library and context
conversion_data_chloroplast <- chloroplast_data %>%
  group_by(sample, context) %>%
  summarise(
    total_methylated = sum(methylated_count),
    total_unmethylated = sum(unmethylated_count),
    total_cytosines = total_methylated + total_unmethylated,
    conversion_rate = 1 - (total_methylated / total_cytosines)
  ) %>%
  ungroup()

# Calculate conversion rate for each library 
conversion_data_chloroplast_contexts_pooled <- chloroplast_data %>%
  group_by(sample) %>%
  summarise(
    total_methylated = sum(methylated_count),
    total_unmethylated = sum(unmethylated_count),
    total_cytosines = total_methylated + total_unmethylated,
    conversion_rate = 1 - (total_methylated / total_cytosines)
  ) %>%
  ungroup()

# View the results
View(conversion_data_chloroplast)
View(conversion_data_chloroplast_contexts_pooled)

# Mitochondrion -----

# Load the chloroplast data
mitochondrion_data <- read_tsv(
  "/path/to/project/WGBS-selected/data/4_methylation/mitochondrion_data.tsv",
  show_col_types = FALSE,
  col_names = c("chromosome", "start", "end", "context", "methyation_level", "methylated_count", 
                "unmethylated_count", "sample")
)

# Calculate conversion rate for each library and context
conversion_data_mitochondrion <- mitochondrion_data %>%
  group_by(sample, context) %>%
  summarise(
    total_methylated = sum(methylated_count),
    total_unmethylated = sum(unmethylated_count),
    total_cytosines = total_methylated + total_unmethylated,
    conversion_rate = 1 - (total_methylated / total_cytosines)
  ) %>%
  ungroup()

# Calculate conversion rate for each library 
conversion_data_mitochondrion_contexts_pooled <- mitochondrion_data %>%
  group_by(sample) %>%
  summarise(
    total_methylated = sum(methylated_count),
    total_unmethylated = sum(unmethylated_count),
    total_cytosines = total_methylated + total_unmethylated,
    conversion_rate = 1 - (total_methylated / total_cytosines)
  ) %>%
  ungroup()

# View the results
View(conversion_data_mitochondrion)
View(conversion_data_mitochondrion_contexts_pooled)

writexl::write_xlsx(list(chloroplast_by_context = conversion_data_chloroplast,
                         chloroplast_context_pooled = conversion_data_chloroplast_contexts_pooled,
                         mitochondrion_by_context = conversion_data_mitochondrion,
                         mitochondrion_context_pooled = conversion_data_mitochondrion_contexts_pooled), 
                    path = "/path/to/project/WGBS-selected/data/4_methylation/bisulfite-converion-rate-selected.xlsx")
