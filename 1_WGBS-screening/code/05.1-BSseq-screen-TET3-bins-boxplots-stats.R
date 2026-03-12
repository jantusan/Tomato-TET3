# Libraries
library(tidyverse)
library(cowplot); theme_set(theme_cowplot())
library(patchwork)
library(ggbeeswarm)
library(scales)

bins_long_table <- read_tsv("../data/4_methylation/tomato_tet3_screening.mC_bins.100kbp.tsv", col_types = cols())  %>% 
  # remove some ANSI carriage movement character
  mutate(chr = str_remove_all(chr, "\\\033\\[K")) %>% 
  # filter controls and samples with almost no reads
  filter( ! sample %in% c("A1_2", "B2_13", "neg","PBS", "Undetermined")) %>% 
  # order the samples as a factor
  mutate(sample = factor(sample, levels = c(
    "WT_1", "WT_2", "WT_3", "WT_4", "WT_6",
    "A_1", "A_4", "A_12", "A_14", "A_15", "A_17", "A_21", "A_22",  
    "A1_4", "A1_6", "A1_7", "A1_8", 
    "A11_3", "A11_5", "A11_8", "A11_10", "A11_11", 
    "B_2", "B_3", "B_4", "B_5", "B_6", "B_8", "B_10", "B_11", 
    "B1_2", "B1_3", "B1_6", "B1_7", "B1_8", "B1_9",
    "B2_1", "B2_5", "B2_9", "B2_10", "B2_14", "B2_16", "B2_18", "B2_20"
  ))) %>% 
  # remove ch09 with the anomaly and organelle genomes
  filter(!(chr %in% c("SL4.0ch09", "Mt", "Pt")))

sample_table <- read_csv("../data/tet3-bs-screen-sample-table.csv", col_types = cols()) %>% 
  mutate(sample = fct_inorder(sample),
         line = fct_inorder(line),
         condition = fct_inorder(condition),
         description = fct_inorder(description))

bins_long_table_filtered <- bins_long_table %>%
  left_join(sample_table, by = "sample") %>% 
  filter(cytosines > 14e6) # c("B2_18", "B2_20", "B_6", "B_11", "WT_6", "B_5", "A_1", "A_12", "B1_2", "B2_14"))

bins_long_table_filtered %>% 
  mutate(coverage = reads_m + reads_n) %>% 
  ggplot(aes(y = coverage, x = sample, fill = condition)) + 
  facet_grid(cols = vars(description), rows = vars(context), scales = "free", space = "free_x") +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.2, point_colour = NA) + 
  geom_boxplot(width = .15, outlier.shape = NA, color = 1) +
  gghalves::geom_half_point(side = "l", shape = 16, size = 0.5, alpha = 0.3, range_scale = .4) +
  scale_fill_manual(values = subline_colors) + 
  scale_y_continuous(trans = "pseudo_log", breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000),
                     labels = label_number(scale_cut = cut_short_scale())) +
  labs(title = "Cytosine coverage for 100 kbp bins", 
       subtitle = "There is a lot of bins with less than 10 Cs, maybe I should filter them out") +
  theme(axis.line=element_blank(), strip.background = element_rect(colour = NA, fill = "lightpink"),
        panel.border = element_rect(colour="black", size=1, fill=NA), legend.position = "none",
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5), axis.title.x = element_blank())

# stats
stats_tbl <- bins_long_table_filtered %>% 
  group_by(line, condition, description, sample, context, genotype) %>% 
  summarise(mR_avg = mean(mR), .groups = "drop")

cg_aov <- stats_tbl %>% 
  filter(context == "CG") %>% 
  aov(mR_avg ~ description, data = .)
summary(cg_aov)

library(multcomp)
library(multcompView)

cg_tukey <- glht(cg_aov, linfct = mcp(description = "Tukey"))
summary(cg_tukey)

cg_letters <- cld(cg_tukey)
cg_letters

get_letters <- function(df){
  fit <- aov(mR_avg ~ description, data = df)
  tuk <- glht(fit, linfct = mcp(description = "Tukey"))
  let <- cld(tuk)$mcletters$Letters
  tibble(description = names(let), letters = unname(let), p_anova = summary(fit)[[1]][["Pr(>F)"]][1])
}

letters_tbl <- stats_tbl %>% 
  group_by(context) %>% 
  group_modify(~ get_letters(.x)) %>% 
  ungroup()

letters_tbl

# By crhomatin state:

hchrom_H3K27ac_100kbp_bins <- read_tsv("/path/to/genome/SL4.0_heterochromatin_100kbp_bins_H3K27ac_chXX.tsv", show_col_types = FALSE) %>% 
  mutate(seqnames = str_c("SL4.0", seqnames),
         start = start - 1) %>% 
  dplyr::rename(chr = seqnames)

bins_long_table_filtered %>% 
  left_join(hchrom_H3K27ac_100kbp_bins, by = c("chr", "start", "end")) %>% 
  filter(!is.na(chromatin)) %>% 
  group_by(line, condition, description, sample, context, genotype, chromatin) %>% 
  summarise(mR_avg = mean(mR), .groups = "drop") %>% 
  ggplot(aes(description, mR_avg, fill = condition)) +
  facet_grid(cols = vars(line), rows = vars(context, chromatin), scales = "free", space = "free_x") +
  geom_boxplot(outlier.colour = NA, color = "black", alpha = 0.3, width = 0.9) +
  geom_quasirandom(aes(shape = genotype), size = 3, dodge.width = 0.75) +
  scale_shape_manual(values = genotype_shapes) +
  scale_fill_manual(values = subline_colors) + 
  guides(fill = "none", shape = guide_legend(override.aes = list(fill = "darkseagreen2", size = 5))) +
  ylab("% mC") +
  theme(legend.position = "bottom", 
        strip.background = element_rect(colour = NA, fill = "lightpink"),
        panel.border = element_rect(colour = "black", size = 1), axis.line = element_blank(),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5), axis.title.x = element_blank()) + 
  labs(title = "Methylation levels of TET3 plants",
       subtitle = "using average for each plant, split by Eu/Heterochromatin", 
       caption = "each dot represents average methylation rate for the whole genome")

# stats by chromatin state
stats_tbl_chr <- bins_long_table_filtered %>% 
  left_join(hchrom_H3K27ac_100kbp_bins, by = c("chr", "start", "end")) %>% 
  filter(!is.na(chromatin)) %>% 
  group_by(line, condition, description, sample, context, genotype, chromatin) %>% 
  summarise(mR_avg = mean(mR), .groups = "drop")

get_letters_chr <- function(df){
  fit <- aov(mR_avg ~ description, data = df)
  tuk <- glht(fit, linfct = mcp(description = "Tukey"))
  let <- cld(tuk)$mcletters$Letters
  tibble(description = names(let), letters = unname(let), p_anova = summary(fit)[[1]][["Pr(>F)"]][1])
}

letters_tbl_chr <- stats_tbl_chr %>% 
  group_by(context, chromatin) %>% 
  group_modify(~ get_letters_chr(.x)) %>% 
  ungroup()

print(letters_tbl_chr, n = 200)
