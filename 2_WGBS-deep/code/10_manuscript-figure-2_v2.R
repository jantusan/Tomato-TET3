
setwd("/path/to/project/WGBS-selected/analysis/")

# Figure 2A Lineplot ch02 ----
## libraries ----
library(tidyverse)
library(cowplot); theme_set(theme_cowplot())
library(patchwork)
library(slider)
library(ggforce)
library(scales)

## functions ----
## Function to plot chr as facets in ggplot
facet_chr <- function(chr = chr, rows = NULL, scales = "free_x", space = "free_x", breaks = 40e6, switch = NULL) {
  list(facet_grid(rows = vars({{ rows }}), cols = vars({{ chr }}), scales = scales, space = space, switch = switch),
       scale_x_continuous(expand = c(0, 0), limits = c(0, NA), breaks = scales::breaks_width(breaks), 
                          labels = label_number(scale_cut = cut_short_scale())),
       xlab("Genomic coordinates (bp)"))
}

## Data and Function to create heterochromatin ribbon
heterochromatin_H3K27ac = read_tsv("/path/to/genome/SL4.0_heterochromatin_H3K27ac_chXX.bed", col_names = c("chr","start","end"), col_types = cols()) %>% 
  mutate(chr = str_c("SL4.0", chr))

SL4_chrs <- str_c("SL4.0ch", str_pad(c(1:8, 10:12), 2, pad = "0"))

geom_heterochromatin <- function(hchrom = heterochromatin_H3K27ac, chr_col = chr, chrs = SL4_chrs) {
  hchrom <- dplyr::rename(hchrom, {{ chr_col }} := chr)
  list(geom_rect(data = hchrom %>% filter(chr %in% chrs),
                 aes(x = NULL, y = NULL, xmin = start, xmax = end, ymin = -Inf, ymax = Inf), 
                 color = NA, fill = "black", alpha = 0.1))
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

# plot function meth
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

## load data ----
bins_long_table_merged <- read_tsv("../data/4_methylation/tomato_tet3_selected.mC_bins.100kbp.merged.tsv", col_types = cols()) 

bins_wide_table_merged <- bins_long_table_merged %>% 
  pivot_wider(id_cols = chr:context, names_from = sample, values_from = mR, names_sort = TRUE)

bins_long_table_sliding_2Mbp <- bins_wide_table_merged %>%
  mutate(WT_mean = WT) %>%
  pivot_longer_diff_to_reference() %>% 
  compute_sliding_stats(values = mR, window = 20) %>% 
  # order the samples as a factor and add description
  mutate(sample = factor(sample, levels = c("WT", "A4", "A17")),
         description = case_when(sample == "WT" ~ "WT", 
                                 sample == "A4" ~ "TET3", 
                                 sample == "A17" ~ "NTS") %>% 
           fct_inorder())

## plot ----
description_colors <- c(WT  = "grey40", "TET3" = "tomato2", "NTS" = "slateblue2")

bins_long_table_sliding_2Mbp  %>% 
  filter(sample != "WT") %>% 
  plot_genome(y = mean_diff, chrs = "SL4.0ch02", points = TRUE, y_points = diff, breaks = 10e6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") +
  geom_line(aes(group = sample), linewidth = 0.7) +
  labs(title = "Methylation change for 100 kbp bins",
       y = "% mC difference to WT ") +
  guides(color = guide_legend(keyheight = 1)) +
  theme(legend.position = c(0.73, 0.41), legend.text.align = 0.5, legend.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1), axis.title.x = element_text(angle = 0)) 

ggsave2("../figures/genome_plots/genome-plot-ch02-diff-Figure-2A.pdf", width = 5, height = 6.5)

# Figure 2B DMR colplot ----
## libraries ----
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(scales)
library(patchwork)

## load data ----
tb_dmrs_filtered <- read_tsv("../data/6_dmrs-2025/table-dmrs-whole-genome_A4_A17_CHG-bins-combined_filtered.tsv")

DMR_numbers <- tb_dmrs_filtered %>%
  mutate(comparison = fct_inorder(comparison)) %>% 
  group_by(comparison, context,  regionType) %>%
  summarise(DMRs = n(), length = sum(width), .groups = "drop") %>% 
  dplyr::rename(kind = regionType) %>%
  mutate(kind = factor(kind, levels = c("gain","loss")),
         comparison = fct_inorder(comparison))

comparison_colors  <- c("NTS"  = "slateblue2", "TET3"  = "tomato2", WT  = "grey40")

## plot ----
ggplot(DMR_numbers %>% filter(kind == "loss"), 
       aes(x = length, y = comparison, fill = comparison, label = str_c(round(length/1000, 0), " kb"))) +
  facet_grid(rows = vars(context), cols = vars(kind), switch = "y") +
  geom_col() +
  #geom_shadowtext(position = position_stack(vjust = 0.5), color = "white") +
  geom_text(position = position_stack(vjust = 0.5)) +
  scale_y_discrete(position = "right") +
  scale_x_reverse(labels = scales::label_number(scale_cut = cut_si("bp")), expand = expansion(mult = c(0.05, 0))) +
  scale_fill_manual(values = comparison_colors) +
  theme(strip.background = element_rect(colour = NA, fill = "papayawhip"),legend.position= "none",
        axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = NA)) +
  
  ggplot(DMR_numbers %>% filter(kind == "gain"), 
         aes(x = length, y = comparison, fill = comparison, label = str_c(round(length/1000, 0), " kb"))) +
  facet_grid(rows = vars(context), cols = vars(kind)) +
  geom_col() +
  #geom_shadowtext(position = position_stack(vjust = 0.5), color = "white") +
  geom_text(position = position_stack(vjust = 0.5)) +
  scale_x_continuous(labels = scales::label_number(scale_cut = cut_si("bp")), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = comparison_colors) +
  theme(strip.background = element_rect(colour = NA, fill = "papayawhip"),legend.position= "none",
        axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = NA)) +
  
  plot_annotation(title = "Differentially Methylated Regions (DMRs)",
                  subtitle = "50 bp Bins (some merged). p-value < 0.01 & number of Cs > 3\nMinimum %mC diff: 60% for CG & CHG; 40% for CHH", 
                  caption = "Labels represent total length of DMRs")

ggsave2("../figures/dmrs-2026-CHG-bins/colplot-DMRs-whole-genome-Figure-2B.pdf", width = 8, height = 5)



# Figure 2C Volcano ----
## libraries ----
library(tidyverse)
library(scales)
library(cowplot); theme_set(theme_cowplot())

## functions ----
# function to display y-axis values
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-1000, Inf))
}

# function to plot volcano plots of DMRs
volcano_dmrs <- function(tb_dmrs, v_lines = vlines, pval_limit = 1e-30) {
  ggplot(tb_dmrs %>% filter(cytosinesCount >= 4), aes(x = (proportion2 - proportion1), y = pValue)) +
    facet_grid(cols = vars(comparison), rows = vars(context)) +
    geom_vline(data = v_lines, aes(xintercept = threshold),  color = "dodgerblue3") +
    geom_vline(data = v_lines, aes(xintercept = -threshold), color = "dodgerblue3") +
    geom_hline(yintercept = 0.01, color = "dodgerblue3") +
    geom_point(shape = 16, size = 0.3, alpha = 0.3) +
    scale_x_continuous(expand = c(0.01, 0.01), breaks = breaks_width(0.4), labels = label_percent(), name = "%mC - WT %mC") +
    scale_y_continuous(trans = reverselog_trans(10), limits = c(1, pval_limit),
                       breaks = breaks_log(n = 6), expand = c(0.01, 0.01), 
                       labels = trans_format("log10", math_format(10^.x))) +
    labs(title = "Volcano plot of DMRs") +
    theme(strip.background = element_rect(colour = NA, fill = "papayawhip"),
          panel.border = element_rect(colour = "black", size = 1), axis.line = element_blank())
}

## load data ----
tb_dmrs_whole_genome <- read_tsv("../data/6_dmrs-2025/table-dmrs-whole-genome_A4_A17_CHG-bins-combined_filtered.tsv", show_col_types = FALSE) %>% 
  mutate(comparison = factor(comparison, levels = c("TET3", "NTS")))

# mpd thresholds
vlines <- tibble(context = c("CG", "CHG", "CHH"),
                 threshold = c(0.5, 0.5, 0.3))

## plot ----
volcano_dmrs(tb_dmrs_filtered, pval_limit = 1e-30)

ggsave(filename = "../figures/dmrs-2026-CHG-bins/volcano-DMRs-zoom-30-Figure-2C.pdf", width = 6, height = 7)

# Figure 2D heatmaps ----
## libraries ----
library(tidyverse)
library(ComplexHeatmap)
ht_opt$message = FALSE
library(tidyHeatmap)
library(patchwork)
library(cowplot); theme_set(theme_cowplot())
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

# Function that takes a ComplexHeatmap list and makes it work with Patchwork syntax
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
description_colors <- c(WT  = "grey40", "NTS" = "slateblue2", "TET3" = "tomato2") %>% fct_inorder() # they have to be in this order, for some reason

tb_dmrs_meth_filtered %>% 
  filter_dmrs_and_pivot(comp = "a4", ctx = "CG") %>% 
  heatmap_dmrs_split_no_scaled(title = "A.T1 TET3 CG\n", show_legends = FALSE) +
  
  tb_dmrs_meth_filtered %>% 
  filter_dmrs_and_pivot(comp = "a17", ctx = "CG") %>% 
  heatmap_dmrs_split_no_scaled(title = "A.T1 NT CG\n") +
  
  tb_dmrs_meth_filtered %>% 
  filter_dmrs_and_pivot(comp = "a4", ctx = "CHH") %>% 
  heatmap_dmrs_split_no_scaled(title = "A.T1 TET3 CHH\n", show_legends = FALSE) +
  
  tb_dmrs_meth_filtered %>% 
  filter_dmrs_and_pivot(comp = "a17", ctx = "CHH") %>% 
  heatmap_dmrs_split_no_scaled(title = "A.T1 NT CHH\n") +
  
  plot_layout(nrow = 2, widths = c(8.5, 11))

ggsave2("../figures/dmrs/heatmaps-DMRs-CG-CHH_Figure-2D.pdf", width = 6, height = 10)

# Figure XN whatever ----
## libraries ----
## functions ----
## load data ----
## plot ----