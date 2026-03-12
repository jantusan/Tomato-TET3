
setwd("/path/to/project/WGBS-selected/analysis/")

# Figure 2 supplement 1 --------------------------------------------------------
# Figure 2 Supplement 1A genome plot ----
## libraries ----
library(tidyverse)
library(cowplot); theme_set(theme_cowplot())
library(patchwork)
library(slider)
library(ggforce)
library(scales)

## functions ----
# Function to plot chr as facets in ggplot
facet_chr <- function(chr = chr, rows = NULL, scales = "free_x", space = "free_x", breaks = 40e6, switch = NULL) {
  list(facet_grid(rows = vars({{ rows }}), cols = vars({{ chr }}), scales = scales, space = space, switch = switch),
       scale_x_continuous(expand = c(0, 0), limits = c(0, NA), breaks = scales::breaks_width(breaks), 
                          labels = label_number(scale_cut = cut_short_scale())),
       xlab("Genomic coordinates (bp)"))
}

# Function to pivot a wide methylation table while comparing to a reference
pivot_longer_diff_to_reference <- function(bins_wide, reference = WT_mean, values = mR) {
  bins_wide %>%
    dplyr::rename(ref = {{ reference }}) %>%
    pivot_longer(-c(chr:context,ref), names_to = "sample", values_to = as.character(ensym(values)),
                 names_transform = list(sample = forcats::fct_inorder)) %>%
    mutate(diff = {{ values }} - ref)
}

# Function to compute sliding function at the middle point of a window
slide_dbl_midwindow <- function(.x, .f, window = 100, ...) {
  slide_dbl(.x, ~.f(.x, na.rm = T), ...,  
            .before = ceiling(window / 2) - 1, # to account for the middle point 
            .after = floor(window / 2))        # otherwise the real window size would be window + 1
}

# for some reasons I cannot combine the above function with the embrace operator {{, so I create these:
slide_dbl_midwindow_mean <- function(.x, window = 100, ...) {
  slide_dbl_midwindow(.x, .f = mean, window = window, ...)
}

slide_dbl_midwindow_median <- function(.x, window = 100, ...) {
  slide_dbl_midwindow(.x, .f = median, window = window, ...)
}

# Function to compute moving averages (updated, window split before and after)
compute_sliding_stats <- function(bins_long_with_diff, values = mR, window = 100) {
  bins_long_with_diff %>% 
    arrange(sample, chr, context, start) %>%
    group_by(sample, chr, context) %>%
    mutate(mean        = slide_dbl_midwindow_mean(.x    = {{ values }}, window = window),
           median      = slide_dbl_midwindow_median(.x = {{ values }}, window = window),
           mean_diff   = slide_dbl_midwindow_mean(.x   = diff,         window = window),
           median_diff = slide_dbl_midwindow_median(.x = diff,         window = window)) %>% 
    ungroup()
}

## Data and Function to create heterochromatin ribbon
heterochromatin_H3K27ac = read_tsv("/path/to/genome/SL4.0_heterochromatin_H3K27ac_chXX.bed", col_names = c("chr","start","end"), col_types = cols())

SL4_chrs <- str_c("SL4.0ch", str_pad(c(1:8, 10:12), 2, pad = "0")) %>% 
  str_remove("SL4.0")

# It works now (it didn't work for some time, I'm not sure why)
geom_heterochromatin <- function(hchrom = heterochromatin_H3K27ac, chr_col = chr, chrs = SL4_chrs) {
  hchrom <- hchrom %>% 
    filter(chr %in% chrs) %>% 
    dplyr::rename({{ chr_col }} := chr)
  list(geom_rect(data = hchrom, aes(x = NULL, y = NULL, xmin = start, xmax = end, ymin = -Inf, ymax = Inf), 
                 color = NA, fill = "black", alpha = 0.1))
}

## load data ----
bins_long_table_merged <- read_tsv("../data/4_methylation/tomato_tet3_selected.mC_bins.100kbp.merged.tsv", col_types = cols())

bins_wide_table_merged <- bins_long_table_merged %>% 
  pivot_wider(id_cols = chr:context, names_from = sample, values_from = mR, names_sort = TRUE)

# Compute the sliding binned data for 20 bins of 100kpb, or 2 Mbp
bins_long_table_sliding_2Mbp <- bins_wide_table_merged %>%
  mutate(WT_mean = WT) %>%
  pivot_longer_diff_to_reference() %>% 
  compute_sliding_stats(values = mR, window = 20) %>% 
  # order the samples as a factor and add description
  mutate(sample = factor(sample, levels = c("WT", "A4", "A17")),
         description = case_when(sample == "WT" ~ "WT", 
                                 sample == "A4" ~ "TET3", 
                                 sample == "A17" ~ "NTS") %>% 
           factor(levels = c("WT", "TET3", "NTS")))

## plot ----
plot_genome <- function(data, y, chr_col = chr, scales = "free", chrs = SL4_chrs, breaks = 20e6, points = FALSE, y_points = NULL) {
  data %>%
    filter({{ chr_col }} %in% chrs) %>%
    ggplot(aes(x = start, y = {{ y }}, color = description)) +
    facet_chr(chr = {{ chr_col }}, rows = context, scales = scales, space = "free_x", breaks = breaks) +
    geom_heterochromatin(chrs = chrs) +
    scale_y_continuous(expand = expansion(mult = 0.05), limits = c(NA, NA)) +
    { if(points) geom_point(aes(y = {{ y_points }}), shape = 16, size = 1, alpha = 0.3) } +
    geom_line(aes(group = sample)) +
    scale_color_manual("Sample", values = description_colors)  +
    theme_cowplot() +
    theme(panel.spacing.x=unit(0.5, "lines"),
          strip.background.y = element_rect(colour = NA, fill = "papayawhip"), strip.background.x = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 1), axis.line = element_blank(),
          axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
          legend.position = "bottom", axis.title.x = element_blank())
}

description_colors <- c(WT  = "grey40", "TET3" = "tomato2", "NTS" = "slateblue2")

bins_long_table_sliding_2Mbp  %>% 
  mutate(chr = str_remove(chr, "SL4.0")) %>% 
  filter(sample != "WT") %>% 
  plot_genome(y = mean_diff) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") +
  labs(title = "DNA Methylation across Tomato chromosomes",
       subtitle = "Sliding window of 2Mbp. Shaded region represents Heterochromatin",
       y = "% mC difference to WT (mean)")

ggsave2("../figures/genome_plots/genome-plot-diff-for-figure-supp-.pdf", width = 16, height = 8)



# Figure 2 Supplement 1B boxplots chromatin ----
## libraries ----
library(tidyverse)
library(cowplot); theme_set(theme_cowplot())
library(patchwork)
library(ggbeeswarm)
library(scales)

## load data ----
sample_table <- read_tsv("../data/sample_table.tsv", col_types = cols()) %>% 
  mutate(sample = fct_inorder(sample),
         condition = fct_inorder(condition),
         replica = as_factor(replica),
         description = fct_inorder(description))

bins_long_table_merged <- read_tsv("../data/4_methylation/tomato_tet3_selected.mC_bins.100kbp.merged.tsv") %>% 
  mutate(description = case_when(sample == "WT" ~ "WT", 
                                 sample == "A4" ~ "TET3", 
                                 sample == "A17" ~ "NTS") %>% 
           fct_inorder())

pivot_longer_diff_to_reference <- function(bins_wide, reference = WT_mean, values = mR, remove_reference = TRUE) {
  bins_wide %>%
    rename(ref = {{ reference }}) %>%
    separate_wider_delim(ref, delim = "fake_dlim", names = c("ref"), cols_remove = remove_reference) %>% 
    pivot_longer(-c(chr:context,ref), names_to = "sample", values_to = as.character(ensym(values)),
                 names_transform = list(sample = forcats::fct_inorder)) %>%
    mutate(diff = {{ values }} - ref)
}

bins_long_table_merged_diff <- bins_long_table_merged %>% 
  pivot_wider(id_cols = chr:context, names_from = sample, values_from = mR, names_sort = TRUE) %>% 
  pivot_longer_diff_to_reference(reference = WT, remove_reference = FALSE) %>% 
  # add description
  mutate(description = case_when(sample == "A4" ~ "TET3", 
                                 sample == "A17" ~ "NTS") %>% 
           fct_inorder())

hchrom_H3K27ac_100kbp_bins <- read_tsv("/path/to/genome/SL4.0_heterochromatin_100kbp_bins_H3K27ac_chXX.tsv", show_col_types = FALSE) %>% 
  mutate(seqnames = str_c("SL4.0", seqnames),
         start = start - 1) %>% 
  dplyr::rename(chr = seqnames)

## plot ----
sample_colors <- c(WT  = "grey40", TET3  = "tomato2", NTS = "slateblue2")

((bins_long_table_merged %>% 
    left_join(hchrom_H3K27ac_100kbp_bins, by = c("chr", "start", "end")) %>% 
    filter(!is.na(chromatin)) %>% 
    ggplot(aes(description, mR, fill = description)) +
    ylab("% mC")) /
    (bins_long_table_merged_diff %>% 
       left_join(hchrom_H3K27ac_100kbp_bins, by = c("chr", "start", "end")) %>% 
       filter(!is.na(chromatin)) %>% 
       ggplot(aes(description, diff, fill = description)) +
       geom_hline(yintercept = 0, linetype = "dashed") +
       stat_summary(geom = "col", fun = mean, color = NA, alpha = 0.9, width = 0.15) +
       coord_cartesian(ylim = c(-10, 1)) + 
       ylab("% mC difference to WT")) &
    facet_grid(cols = vars(context, chromatin)) &
    geom_violin(scale = "count", alpha = 0.5, color = NA) &
    geom_boxplot(outlier.shape = NA, coef = 0, notch = TRUE, color = "black", alpha = 0.3, width = 0.3) &
    stat_summary(geom = "point", fun = mean, shape = 23, color = "black", size = 3) &
    scale_fill_manual(values = sample_colors) &
    theme(legend.position = "none", strip.background = element_rect(colour = NA, fill = "lightpink"),
          panel.border = element_rect(colour = "black", linewidth = 1), axis.line = element_blank(),
          axis.title.x = element_blank())) +
  plot_annotation(title = "Methylation levels of TET3 plants by Chromatin compartment",
                  subtitle = "using averages over 100 kbp bins",
                  caption = "Dots represent average of methylation rate")

ggsave2("../figures/mC_bins/boxplot-merged-100kb-chromatin_for-figure-supplement.pdf", width = 11, height = 7)


# Figure 2 supplement 2 --------------------------------------------------------
# Figure 2 supplement 2A coverage ----
## libraries ----
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(scales)
library(patchwork)
library(slider)

## load data ----
cytosine_coverage_tbl_pooled <- read_tsv("../data/4_methylation/ch06/cytosine_coverage_ch06_tbl_pooled.tsv", col_types = cols()) %>% 
  mutate(sample = toupper(sample),
         sample = case_when(sample == "WT" ~ "WT", 
                            sample == "A4" ~ "TET3", 
                            sample == "A17" ~ "NTS") %>% 
           fct_inorder())

## plot ----
sample_colors <- c(WT  = "grey40", TET3  = "tomato2", NTS = "slateblue2")

coverage_base_plot_pooled <- cytosine_coverage_tbl_pooled %>% 
  group_by(sample, context) %>% 
  mutate(total = sum(count),
         freq = count / total, 
         left = slide_dbl(freq, ~ 1 - sum(.), .before = Inf)) %>% 
  ggplot(aes(x = readsN, color = sample, fill = sample)) +
  facet_wrap(vars(context)) +
  scale_color_manual(values = sample_colors) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(limits = c(0, NA), labels = scales::label_percent()) +
  coord_cartesian(xlim = c(0, 20)) +
  xlab("read coverage") +
  theme(strip.background = element_rect(fill = "darkseagreen1"),
        panel.border = element_rect(colour = "black", linewidth = 1), axis.line = element_blank())

(coverage_base_plot_pooled +
    geom_line(aes(y = freq), linewidth = 2, alpha = 0.5) +
    geom_point(aes(y = freq), shape = 21, size = 2, color = 1) +
    ylab("frequency"))  /
  (coverage_base_plot_pooled +
     geom_line(aes(x = readsN, y = left), linewidth = 2, alpha = 0.5) +
     geom_point(aes(x = readsN, y = left), shape = 21, size = 2, color = 1) +
     ylab("right tail")) +
  plot_annotation(title = "Read coverage per cytosine", 
                  subtitle = "for chromosome 6 of tomato (47 Mbp)\nReplicas pooled",
                  caption = "Top: proportion of the cytosines with a certain coverage
                             Bottom: proportion of the cytosines over a certain coverage threshold")  +
  plot_layout(guides = "collect")

ggsave(filename = "../figures/ch06-test/meth-stats/cytosine-coverage-ch06-freq-pooled_for-figure-supplement.pdf", height = 7, width = 7)

# Figure 2 supplement 2B correlation ----
## libraries ----
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(scales)
library(patchwork)
library(slider)

## functions ----
## load data ----
spatial_correlation_pooled <- read_tsv(file = "../data/4_methylation/ch06/spatial_correlation_ch06_pooled.tsv", col_types = cols()) %>% 
  mutate(sample = toupper(sample),
         sample = case_when(sample == "WT" ~ "WT", 
                            sample == "A4" ~ "TET3", 
                            sample == "A17" ~ "NTS") %>% 
           fct_inorder())

## plot ----
spatial_correlation_plot_linetype_pooled <- spatial_correlation_pooled %>% 
  ggplot(aes(x = distance, y = correlation, color = context, fill = context, group = str_c(context, sample), shape = sample)) +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
  geom_line(linewidth = 2, alpha = 0.5) +
  geom_point(size = 2, color = 1) +
  xlab("Distance (bp)") +
  theme(strip.background = element_blank(), panel.border = element_rect(colour = "black", linewidth = 1),
        axis.line = element_blank())

(spatial_correlation_plot_linetype_pooled + coord_cartesian(xlim = c(0, 500)) + theme(axis.text.x = element_text(angle = -45, hjust = 0))) +
  spatial_correlation_plot_linetype_pooled + 
  scale_x_continuous(trans='log10', breaks = scales::breaks_log(), labels = label_number(scale_cut = cut_si(""))) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Spatial correlation of DNA methylation values",
                  subtitle = "for chromosome 6 of tomato (47 Mbp)\nReplicas pooled",
                  caption = "left: linear scale, right: semilogarithmic scale")

ggsave(filename = "../figures/ch06-test/meth-stats/spatial-correlation-ch06-pooled_for-figure-supplement.pdf", height = 5, width = 9)

# Figure 2 supplement 2 --------------------------------------------------------
# Figure 2 supplement 2C DMRs bins ----
## libraries ----
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(scales)

## load data ----
stats_dmrs_bins <- read_tsv("../data/4_methylation/ch06/dmrs_ch06_a17_a4_bins_size-stats.tsv", show_col_types = FALSE)  %>% 
  mutate(sample = toupper(sample),
         sample = case_when(sample == "WT" ~ "WT", 
                            sample == "A4" ~ "TET3", 
                            sample == "A17" ~ "NTS") %>% 
           factor(levels = c("WT", "TET3", "NTS")))

## plot ----
ggplot(stats_dmrs_bins, aes(x = as_factor(size), y = as_factor(min_prop), fill = length)) +
  facet_grid(cols = vars(sample, regionType), rows = vars(context), space = "free") +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#340240") +
  geom_tile() +
  scale_fill_viridis_c("Length\n(bp)", trans='log10', breaks = scales::breaks_log(), labels = label_number(scale_cut = cut_short_scale())) +
  #scale_fill_viridis_c() +
  labs(title = "DMRs with Bins method",
       subtitle = "for chromosome 6 of tomato (47 Mbp)", 
       x = "Moving window size",
       y = "Minimun proportion difference") +
  theme(legend.key.height = unit(0.6, "cm"), strip.background = element_rect(colour = NA, fill = "lavender"),
        panel.border = element_rect(colour = "black", linewidth = 1), axis.line = element_blank(),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) 

ggsave2("../figures/ch06-test/dmrs-bins/dmrs-bins-heatmap-length_for-figure-supplement.pdf", height = 6.5, width = 8)

# Figure 2 supplement 2C DMRs noise ----
## libraries ----
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(scales)

## load data ----
stats_dmrs_noise <- read_tsv("../data/4_methylation/ch06/dmrs_ch06_a17_a4_noise_windows-stats.tsv", show_col_types = FALSE)  %>% 
  mutate(sample = toupper(sample),
         sample = case_when(sample == "WT" ~ "WT", 
                            sample == "A4" ~ "TET3", 
                            sample == "A17" ~ "NTS") %>% 
           factor(levels = c("WT", "TET3", "NTS")))

## plot ----
ggplot(stats_dmrs_noise, aes(x = as_factor(size), y = as_factor(min_prop), fill = length)) +
  facet_grid(cols = vars(sample, regionType), rows = vars(context), space = "free") +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#340240") +
  geom_tile() +
  scale_fill_viridis_c("Length\n(bp)", trans='log10', breaks = scales::breaks_log(), labels = label_number(scale_cut = cut_short_scale())) +
  #scale_fill_viridis_c() +
  labs(title = "DMRs with Noise filtering method",
       subtitle = "for chromosome 6 of tomato (47 Mbp)", 
       x = "Moving window size",
       y = "Minimun proportion difference") +
  theme(legend.key.height = unit(0.6, "cm"), strip.background = element_rect(colour = NA, fill = "lavender"),
        panel.border = element_rect(colour = "black", linewidth = 1), axis.line = element_blank(),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) 

ggsave2("../figures/ch06-test/dmrs-noise-filter/dmrs-noise-filter-heatmap-length_for-figure-supplement.pdf", height = 6.5, width = 8)

# Figure 2 supplement 3 CHG heatmap ----
## libraries ----
library(tidyverse)
library(ComplexHeatmap)
ht_opt$message = FALSE
library(tidyHeatmap)
library(patchwork)
library(cowplot): theme_set(theme_cowplot())
library(scales)

## functions ----
# color ramp from a Brewer palette
brewer.ramp <- function(start, end, pal) {
  palette = RColorBrewer::brewer.pal(n = RColorBrewer::brewer.pal.info[pal, "maxcolors"],  name = pal)
  sequence = seq(start, end, (end - start)/(RColorBrewer::brewer.pal.info[pal, "maxcolors"] - 1))
  ramp = circlize::colorRamp2(sequence, palette)
  return(ramp)
}

# Function that takes a tibble of combined DMR comparisons and returns a long table
filter_dmrs_and_pivot <- function(tb_dmrs, comp, ctx) {
  tb_dmrs %>% 
    filter(comparison == comp, 
           context == ctx)  %>% 
    pivot_longer(cols = starts_with("proportion_"), values_to = "methylation", names_to = "sample", 
                 names_prefix = "proportion_", names_transform = toupper) %>% 
    mutate(index = str_c(seqnames, ":", start, ",", end),
           # add description
           description = case_when(sample == "WT" ~ "WT", 
                                   sample == "A4" ~ "TET3", 
                                   sample == "A17" ~ "NTS") %>% 
             fct_inorder()) %>% 
    # scale methylation
    group_by(index) %>% 
    mutate(scaled = scale(methylation)) %>% 
    ungroup()
}

#Function that takes a ComplexHeatmap list and makes it work with Patchwork syntax
wrap_heatmap_list <- function(heatmap_list, draw = FALSE) {
  heatmap_list %>% 
    {if (draw) ComplexHeatmap::draw(.) else .} %>%
    grid::grid.grabExpr() %>% 
    patchwork::wrap_elements()
}

heatmap_dmrs_split_no_scaled <- function(tb_dmrs_meth_long, show_legends = TRUE, title = character(0), max_n = 4000, 
                                         top_tile_variable = description, top_tile_colors = description_colors, subset_trick = FALSE,
                                         use_raster = FALSE) {
  
  n_dmrs = unique(tb_dmrs_meth_long$index) %>% length()
  
  if (subset_trick == TRUE) {
    # remove some DMRs, as if not, the diff tile doesn't match the heatmap
    tb_dmrs_meth_long <- tb_dmrs_meth_long %>% 
      filter(seqnames %in% str_c("chr", 1:9))
  }
  
  final_n_dmrs = unique(tb_dmrs_meth_long$index) %>% length()
  if (final_n_dmrs < max_n) max_n <- final_n_dmrs
  
  tb_dmrs_meth_long_sample <- tb_dmrs_meth_long %>% 
    filter(index %in% sample(unique(tb_dmrs_meth_long$index), size = max_n)) %>% 
    group_by(regionType) %>% 
    mutate(genotype = {{ top_tile_variable }})
  
  (heatmap(tb_dmrs_meth_long_sample, .row = index, .column = description, .value = methylation, 
           cluster_columns = TRUE, show_row_names = FALSE, use_raster = use_raster, 
           palette_value = color_ramp_unscaled, name = "mC", palette_grouping = list(c("firebrick3", "dodgerblue3")), 
           row_title = NULL, column_title = NULL, width = unit(3, "cm")) %>%
      annotation_tile(genotype, palette = top_tile_colors, show_annotation_name = FALSE) %>% 
      annotation_tile(diff, palette = color_ramp_diff, annotation_name_side = "top", size = unit(1, "cm"))) %>% # + 
    #heatmap(tb_dmrs_meth_long_sample, .row = index, .column = description, .value = scaled, 
    #        cluster_columns = TRUE, show_row_names = FALSE, use_raster = use_raster, 
    #        palette_value = color_ramp_scaled, row_km = 2, name = "scaled", palette_grouping = list(c("firebrick3", "dodgerblue3"), use_raster = TRUE), 
    #        row_title = NULL, column_title = NULL) %>% 
    #annotation_tile(genotype, palette = top_tile_colors, show_annotation_name = FALSE)) %>% 
    as_ComplexHeatmap() %>% 
    draw(column_title = str_c(title, " (" , n_dmrs, " DMRs)"), merge_legend = TRUE, 
         show_heatmap_legend = show_legends, show_annotation_legend = show_legends) %>% 
    wrap_heatmap_list(draw = FALSE)
}

## load data ----
tb_dmrs_meth_filtered <- read_tsv("../data/6_dmrs-2025/table-dmrs-whole-genome-a17-a4-CHG-bins-combined-all-samples-filtered.tsv", show_col_types = FALSE)

## plot ----
color_ramp_diff  <- brewer.ramp(1, -1, "Spectral")
color_ramp_scaled   <- brewer.ramp(2, -2, "RdBu")
color_ramp_unscaled <- brewer.ramp(0, 1, "YlOrRd")
description_colors <- c(WT  = "grey40", "NTS" = "slateblue2", "TET3" = "tomato2") %>% fct_inorder()

tb_dmrs_meth_filtered_combined %>% 
  filter_dmrs_and_pivot(comp = "a4", ctx = "CHG") %>% 
  heatmap_dmrs_split_no_scaled(title = "TET3 CHG\n", show_legends = FALSE) +
  
  tb_dmrs_meth_filtered_combined %>% 
  filter_dmrs_and_pivot(comp = "a17", ctx = "CHG") %>% 
  heatmap_dmrs_split_no_scaled(title = "NTS CHG\n") +
  
  plot_layout(nrow = 1, widths = c(8.5, 11))

ggsave2("../figures/dmrs-2026-CHG-bins/heatmaps-DMRs-CHG_Figure-2S3.pdf", width = 6, height = 5)

# Figure X supplement X --------------------------------------------------------
# Figure X supplement XA ----
## libraries ----
## functions ----
## load data ----
## plot ----