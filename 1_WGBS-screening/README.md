# 1_WGBS-screening

This folder contains the WGBS screening analysis stream for the tomato TET3 study.

## Contents

- `code/`: Screening analysis notebooks and supporting scripts (`01` to `09`)
- `reports/`: Rendered `.nb.html` reports corresponding to notebook analyses
- `data/`: Processed screening-level methylation/chromatin tables

## Figure traceability

- Code series `01` to `06` and `09` provides the screening analysis outputs used for manuscript screening-related figures and supplementary panels.
- Supporting statistical calculations for screening boxplots are in:
  - `code/05.1-BSseq-screen-TET3-bins-boxplots-stats.R`

## Final tables

- Core screening tables are stored under `data/` (e.g., 100 kbp methylation bin summaries and sample tables).
