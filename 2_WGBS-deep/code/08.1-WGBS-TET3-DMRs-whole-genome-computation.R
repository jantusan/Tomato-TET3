# Bash work ----
# on Avon. be inside a screen 
cat("
# request an interactive session on the hmem partition, compute partition gets killed
salloc --nodes=1 --ntasks=24 --mem-per-cpu=32000 --partition=hmem --time=05:30:00
# once inside the interactive
module load GCC/13.2.0
module load OpenMPI/4.1.6
module load R/4.3.3
R
")

# Libraries ----
library(DMRcaller)
library(tictoc)
library(tibble)
library(readr)
#library(plyranges)
#library(betareg)
library(dplyr)
library(purrr)

setwd("/path/to/project/data/4_methylation")

# load data ----
tic()
load("meth_all_replicas_pooled.Rdata")
toc()

# CHG bins -----
## Compute DMRs -----

# New relaxed parameters: bins, 100 bp, mpd 0.5
dmr_params <- list(
  method   = c(CHG = "bins"),
  pvalue   = c(CHG = 0.01),
  size     = c(CHG = 100),
  min_prop = c(CHG = 0.5)
)

compute_dmrs <- function(methylation_data_1, methylation_data_2, context, min_prop, method = "bins", 
                         size = 100, pvalue = 0.01, min_size = 50, min_gap = (2 * size)) {
  computeDMRs(methylationData1 = methylation_data_1, methylationData2 = methylation_data_2,
              context = context, method = method, windowSize = size, binSize = size, 
              minSize = min_size, minGap = min_gap, test = "score", pValueThreshold = pvalue, 
              minCytosinesCount = 1, minProportionDifference = min_prop, minReadsPerCytosine = 1, cores = 16)
}


tic.clearlog()
dmrs_whole_genome <- list()
for (sample in c("a4", "a17")) {
  for (ctx in c("CHG")) {
    tic(paste(sample, ctx))
    dmrs_whole_genome[[sample]][[ctx]] <- compute_dmrs(
      methylation_data_1 = CX_pooled$wt, 
      methylation_data_2 = CX_pooled[[sample]],
      context  = ctx, 
      method   = dmr_params$method[[ctx]],
      pvalue   = dmr_params$pvalue[[ctx]],
      size     = dmr_params$size[[ctx]], 
      min_prop = dmr_params$min_prop[[ctx]])
    toc(log = TRUE)
  }
}
tic.log() %>% unlist() %>% as_tibble() %>% write_tsv("dmrs_CHG_relaxed_whole_genome_timing.tsv")

saveRDS(object = dmrs_whole_genome, file = "../6_dmrs/dmrs_whole_genome_A4_A17_CHG_relaxed.rds")

## Extract mC all samples -----
dmrs_whole_genome <- readRDS(file = "../6_dmrs/dmrs_whole_genome_A4_A17_CHG_relaxed.rds")

tic.clearlog()
dmr_status <- list()
for (dmrs in c("a4", "a17")) {
  for (ctx in c("CHG")) {
    dmr_status[[dmrs]][[ctx]] <- dmrs_whole_genome[[dmrs]][[ctx]]
    for (sample in c("wt", "a17", "a4")) {
      tic(paste(dmrs, ctx, sample))
      dmr_status[[dmrs]][[ctx]] = analyseReadsInsideRegionsForCondition(
        regions = dmr_status[[dmrs]][[ctx]], 
        methylationData = CX_pooled[[sample]], 
        context = ctx, 
        label = paste0("_", sample, "_")
      )
      toc(log = TRUE)
    }
  }
}
tic.log() %>% unlist() %>% as_tibble() %>% write_tsv("dmrs-CHG-bins-add-samples-data-timing.tsv")

save(dmr_status, file = "../6_dmrs/dmrs_whole_genome_A4_A17_CHG_relaxed-all-samples-data.Rdata")


# CG bins -----
## Compute DMRs -----

# New relaxed parameters: bins, 100 bp, mpd 0.5
dmr_params <- list(
  method   = c(CG = "bins"),
  pvalue   = c(CG = 0.01),
  size     = c(CG = 100),
  min_prop = c(CG = 0.5)
)

compute_dmrs <- function(methylation_data_1, methylation_data_2, context, min_prop, method = "bins", 
                         size = 100, pvalue = 0.01, min_size = 50, min_gap = (2 * size)) {
  computeDMRs(methylationData1 = methylation_data_1, methylationData2 = methylation_data_2,
              context = context, method = method, windowSize = size, binSize = size, 
              minSize = min_size, minGap = min_gap, test = "score", pValueThreshold = pvalue, 
              minCytosinesCount = 1, minProportionDifference = min_prop, minReadsPerCytosine = 1, cores = 16)
}


tic.clearlog()
dmrs_whole_genome <- list()
for (sample in c("a4", "a17")) {
  for (ctx in c("CG")) {
    tic(paste(sample, ctx))
    dmrs_whole_genome[[sample]][[ctx]] <- compute_dmrs(
      methylation_data_1 = CX_pooled$wt, 
      methylation_data_2 = CX_pooled[[sample]],
      context  = ctx, 
      method   = dmr_params$method[[ctx]],
      pvalue   = dmr_params$pvalue[[ctx]],
      size     = dmr_params$size[[ctx]], 
      min_prop = dmr_params$min_prop[[ctx]])
    toc(log = TRUE)
  }
}
tic.log() %>% unlist() %>% as_tibble() %>% write_tsv("dmrs_CG_relaxed_whole_genome_timing.tsv")

saveRDS(object = dmrs_whole_genome, file = "../6_dmrs/dmrs_whole_genome_A4_A17_CG_relaxed.rds")

## Extract mC all samples -----
dmrs_whole_genome <- readRDS(file = "../6_dmrs/dmrs_whole_genome_A4_A17_CG_relaxed.rds")

tic.clearlog()
dmr_status <- list()
for (dmrs in c("a4", "a17")) {
  for (ctx in c("CG")) {
    dmr_status[[dmrs]][[ctx]] <- dmrs_whole_genome[[dmrs]][[ctx]]
    for (sample in c("wt", "a17", "a4")) {
      tic(paste(dmrs, ctx, sample))
      dmr_status[[dmrs]][[ctx]] = analyseReadsInsideRegionsForCondition(
        regions = dmr_status[[dmrs]][[ctx]], 
        methylationData = CX_pooled[[sample]], 
        context = ctx, 
        label = paste0("_", sample, "_")
      )
      toc(log = TRUE)
    }
  }
}
tic.log() %>% unlist() %>% as_tibble() %>% write_tsv("dmrs-CG-bins-add-samples-data-timing.tsv")

save(dmr_status, file = "../6_dmrs/dmrs_whole_genome_A4_A17_CG_relaxed-all-samples-data.Rdata")
