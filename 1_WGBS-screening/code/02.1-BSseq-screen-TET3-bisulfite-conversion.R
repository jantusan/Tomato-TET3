library(tidyverse)

# Chloroplast -----

# Load the chloroplast data
chloroplast_data <- read_tsv(
  "/path/to/project/WGBS-screening/data/4_methylation/chloroplast_data.tsv",
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
  "/path/to/project/WGBS-screening/data/4_methylation/mitochondrion_data.tsv",
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
                    path = "/path/to/project/WGBS-screening/data/4_methylation/bisulfite-converion-rate-screening.xlsx")

conversion_data_mitochondrion_contexts_pooled %>% 
  filter(sample %in% c(
    "WT_1", "WT_2", "WT_3", "WT_4", "WT_6",
    "A_1", "A_4", "A_12", "A_14", "A_15", "A_17", "A_21", "A_22",  
    "A1_4", "A1_6", "A1_7", "A1_8", 
    "A11_3", "A11_5", "A11_8", "A11_10", "A11_11", 
    "B_2", "B_3", "B_4", "B_5", "B_6", "B_8", "B_10", "B_11", 
    "B1_2", "B1_3", "B1_6", "B1_7", "B1_8", "B1_9",
    "B2_1", "B2_5", "B2_9", "B2_10", "B2_14", "B2_16", "B2_18", "B2_20"
  )) %>% 
  mutate(sample = factor(sample, levels = c(
    "WT_1", "WT_2", "WT_3", "WT_4", "WT_6",
    "A_1", "A_4", "A_12", "A_14", "A_15", "A_17", "A_21", "A_22",  
    "A1_4", "A1_6", "A1_7", "A1_8", 
    "A11_3", "A11_5", "A11_8", "A11_10", "A11_11", 
    "B_2", "B_3", "B_4", "B_5", "B_6", "B_8", "B_10", "B_11", 
    "B1_2", "B1_3", "B1_6", "B1_7", "B1_8", "B1_9",
    "B2_1", "B2_5", "B2_9", "B2_10", "B2_14", "B2_16", "B2_18", "B2_20"
  ))) %>% 
  arrange(sample) %>% 
  print(n = 50)
