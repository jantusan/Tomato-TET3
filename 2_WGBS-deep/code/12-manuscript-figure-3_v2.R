
setwd("/path/to/project/WGBS-selected/analysis/")

# Figure 3A DMRs genome ----
## libraries ----
library(tidyverse)
library(ggforce)
library(cowplot); theme_set(theme_cowplot())
library(scales)
library(ggblend)
library(ggrastr)

## functions ----
# Function to plot chr as facets in ggplot
facet_chr <- function(chr = chr, rows = NULL, scales = "free_x", space = "free_x", breaks = 40e6, switch = NULL) {
  list(facet_grid(rows = vars({{ rows }}), cols = vars({{ chr }}), scales = scales, space = space, switch = switch),
       scale_x_continuous(expand = c(0, 0), limits = c(0, NA), breaks = scales::breaks_width(breaks), 
                          labels = label_number(scale_cut = cut_short_scale())),
       xlab("Genomic coordinates (bp)"))
}

## Data and Function to create heterochromatin ribbon
heterochromatin_H3K27ac = read_tsv("/path/to/genome/SL4.0_heterochromatin_H3K27ac_chXX.bed", col_names = c("chr","start","end"), col_types = cols())

SL4_chrs <- str_c("ch", str_pad(c(1:8, 10:12), 2, pad = "0"))

geom_heterochromatin <- function(hchrom = heterochromatin_H3K27ac, chr_col = chr, chrs = SL4_chrs) {
  hchrom <- dplyr::rename(hchrom, {{ chr_col }} := chr)
  list(geom_rect(data = hchrom %>% filter(chr %in% chrs),
                 aes(x = NULL, y = NULL, xmin = start, xmax = end, ymin = -Inf, ymax = Inf), 
                 inherit.aes = FALSE, color = NA, fill = "black", alpha = 0.1))
}

## load data ----
tb_dmrs_filtered <- read_tsv("../data/6_dmrs-2025/table-dmrs-whole-genome_A4_A17_CHG-bins-combined_filtered.tsv", show_col_types = FALSE) %>% 
  mutate(comparison = factor(comparison, levels = c("TET3", "NTS")))

## plot ----
comparison_colors  <- c("NTS"  = "slateblue2", "TET3"  = "tomato2", WT  = "grey40")

tb_dmrs_filtered %>% 
  filter(chr %in% SL4_chrs) %>% 
  mutate(mC_diff = proportion2 - proportion1) %>% 
  ggplot(aes(x = start, y = mC_diff, color = comparison, partition = comparison)) +
  #facet_chr(chr = seqnames, rows = (context, regionType), scales = "free", space = "free_x", breaks = 20e6) +
  facet_grid(rows = vars(context, regionType), cols = vars(chr), scales = "free", space = "free") +
  geom_heterochromatin() +
  #geom_point(shape = 16, size = 0.5, alpha = 0.6) |> partition(vars(description)) |> blend("multiply") +
  #geom_point(shape = 16, size = 0.6, alpha = 0.4) * (blend("darken") + blend("multiply", alpha = 0.5)) +
  geom_point(shape = 16, size = 0.6, alpha = 0.4) %>% rasterise(dpi = 300) * 
     (blend("darken") + blend("multiply", alpha = 0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA), breaks = scales::breaks_width(20e6), 
                     labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(expand = c(0, 0), breaks = breaks_width(0.2), labels = label_percent(), name = "%mC - WT %mC") +
  scale_color_manual("comparison", values = comparison_colors)  +
  labs(title = "DMR methylation difference across the genome") +
  xlab("Genomic coordinates (bp)") +
  theme(panel.spacing.x=unit(0.5, "lines"),
        strip.background.y = element_rect(colour = NA, fill = "papayawhip"), strip.background.x = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1), axis.line = element_blank(),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
        legend.position = "bottom", axis.title.x = element_blank())

ggsave("../figures/dmrs-2026-CHG-bins/genome-plot-DMRs-diff-blended-Figure-3A.pdf", width = 12, height = 8, dpi = 600, bg = "white", device = cairo_pdf)

# Figure 3B correlation DMRs genes ----
## libraries ----
library(tidyverse)
library(ComplexHeatmap)
library(shadowtext)

## functions ----
brewer.ramp <- function(start, end, pal) {
  palette = RColorBrewer::brewer.pal(n = RColorBrewer::brewer.pal.info[pal, "maxcolors"],  name = pal)
  sequence = seq(start, end, (end - start)/(RColorBrewer::brewer.pal.info[pal, "maxcolors"] - 1))
  ramp = circlize::colorRamp2(sequence, palette)
  return(ramp)
}

## load data ----
tb_dmrs_filtered <- read_tsv("../data/6_dmrs-2025/table-dmrs-whole-genome_A4_A17_CHG-bins-combined_filtered.tsv", show_col_types = FALSE) %>% 
  mutate(comparison = factor(comparison, levels = c("TET3", "NTS")))

genes = read_rds("/path/to/genome/v4_Mt_Pt_TET3/ITAG4.1_gene_models_Mt_Pt_TET3-named.RDS")

## compute correlation ----
bin_width = 2e6

binned_genes <- genes %>%
  as_tibble() %>%
  group_by(seqnames, bins = cut_width(start, width = bin_width, boundary = 0)) %>%
  summarise(n_genes = n(), .groups = "drop") %>%
  mutate(n_genes = 100 * n_genes / length(genes))

binned_dmrs <- tb_dmrs_filtered %>%
  #filter(context != "CHG", regionType == "gain") %>% 
  group_by(context, regionType, comparison) %>%
  mutate(total = n()) %>% 
  group_by(context, seqnames, regionType, comparison, bins = cut_width(start, width = bin_width, boundary = 0), total) %>%
  summarise(dmrs = n(), .groups = "drop") %>%
  mutate(dmrs = 100 * dmrs / total)
  
binned_joined <- full_join(binned_genes, binned_dmrs,  by = c("seqnames", "bins")) %>%
  na.omit() 
  
# Calculate correlation coefficients for each combination
correlation_data <- binned_joined %>%
  group_by(context, regionType, comparison) %>%
  #filter(n() > 10) %>% 
  summarize(cor = cor(n_genes, dmrs, use = "complete.obs"), .groups = "drop")

# Spread the data for heatmap plotting
heatmap_data <- correlation_data %>%
  pivot_wider(names_from = c(regionType, comparison), values_from = cor) %>% 
  column_to_rownames("context")

## plot heatmap ----
mat = as.matrix(heatmap_data)

# make my own heatmap function
hms = Heatmap(mat, col =  brewer.ramp(1, -1, "RdYlBu"), name = "R ", rect_gp = gpar(col = "black", lwd = 2),
              column_title = "R coefficient", cluster_columns = FALSE, cluster_rows = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                #grid.shadowtext(sprintf(sprintf, mat2[i, j]), x, y, gp = gpar(fontsize = 13))
                grid.shadowtext(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 12 , col = "white"), bg.colour = "black")
              }
)

# draw the heatmaps and save them
draw(hms, column_title = gt_render("**Correlations DMRs vs genes**"))
dev.copy2pdf(file = "../figures/dmrs-2026-CHG-bins/correlation-DMRs-genes-heatmap-2-Mbp_Figure-3B.pdf", width = 3.5, height = 3) %>% invisible()

# Figure 3C metagene plot of DMR location ----
## functions ----
source("/path/to/project/bin/R/enriched_heatmaps_functions.R")
library(plyranges)
library(cowplot); theme_set(theme_cowplot())

## load data ----
tb_dmrs_filtered <- read_tsv("../data/6_dmrs-2025/table-dmrs-whole-genome_A4_A17_CHG-bins-combined_filtered.tsv", show_col_types = FALSE) %>% 
  mutate(comparison = factor(comparison, levels = c("TET3", "NTS")))

genes = read_rds("/path/to/genome/v4_Mt_Pt_TET3/ITAG4.1_gene_models_Mt_Pt_TET3-named.RDS")

test <- tb_dmrs_filtered %>% 
  filter(comparison == "NTS", 
         context == "CHG", 
         regionType == "loss") %>% 
  as_granges() %>% 
  filter_by_overlaps(genes, maxgap = 2000)

meta_tb <- tibble()
for (comp in c("TET3", "NTS")) {
  for (ctx in c("CG", "CHG", "CHH")) {
    for (kind in c("gain", "loss")) {
      dmrs <- read_gff3(str_c("../data/6_dmrs-2025/gffs/DMRs_", comp, "_", ctx, "_", kind, ".gff3"))
      
      mat <- normalise_to_matrix(signal = dmrs, target = filter_by_overlaps(genes, dmrs, maxgap = 2000), 
                                 bed = TRUE, target_ratio = 1/3) 
      meta_tb <- mat[] %>% 
        as_tibble() %>% 
        replace(is.na(.), 0) %>%
        pivot_longer(everything(), names_to = "position", values_to = "value") %>%
        mutate(position = factor(position, levels = unique(position))) %>%
        group_by(position) %>%
        summarise(mean   =   mean(value, na.rm = T), .groups = "drop") %>%
        mutate(r_mean   = slide_dbl(mean,   ~mean(.x,   na.rm = TRUE), .before = 2, .after = 2), 
               comparison = comp, 
               context = ctx, 
               kind = kind) %>% 
        bind_rows(meta_tb)
      
      
    }
  }
}

## plot ----
sample_colors <- c(WT  = "grey40", TET3 = "tomato2", NTS = "slateblue2")

ggplot(meta_tb, aes(x = position, y = mean, color = comparison, linetype = kind, group = paste(comparison, kind))) +
  facet_grid(rows = vars(context), scales = "free_y") +
  #geom_point(pch=16, size=1, alpha = 0.6) +
  geom_line(mapping = aes(y = r_mean), linewidth = 0.5, alpha=1) +
  scale_x_discrete(breaks = c("u1","t1","d1","d50"), labels = c(paste0("-",2," Kb"), "Start", "End", paste0("+",2," Kb"))) +
  scale_y_continuous("genes with DMRs", labels = scales::label_percent()) +
  scale_color_manual(values = sample_colors) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1) + # y axis on all facets
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1) + # y axis on all facets
  annotate("rect", xmin="t1", xmax="d1", ymin=-Inf, ymax=Inf, alpha=0.1, fill="black") +
  theme(axis.line=element_blank(), axis.title.x = element_blank(),
        strip.background.y = element_rect(colour = NA, fill = "papayawhip"), strip.background.x = element_blank(),
        panel.spacing=unit(2, "lines"), strip.background=element_blank()) +
  labs(title = "Distribution of DMRs across genes")

ggsave2("../figures/dmrs-2026-CHG-bins/metaplot-DMRs-genes_Figure-3C-.pdf", width = 120, height = 100, units = "mm")
