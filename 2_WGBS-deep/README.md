# 2_WGBS-deep

This folder contains deep WGBS analyses, DMR workflows, and manuscript figure-generation scripts for the deep methylome section.

## Contents

- `code/`: Deep WGBS notebooks (`01` to `09`) and figure/computation scripts
- `reports/`: Rendered `.nb.html` reports for each deep-analysis notebook
- `data/`: Processed DMR/correlation/coverage tables used in notebooks and figure scripts

## Figure traceability

- Manuscript **Figure 2** panels are generated in:
  - `code/10_manuscript-figure-2_v2.R`
- Manuscript **Figure 2 supplements** are generated in:
  - `code/11_manuscript-figure-2-supplements.R`
- Manuscript **Figure 3** panels are generated in:
  - `code/12-manuscript-figure-3_v2.R`
- Whole-genome DMR computation feeding downstream plots is in:
  - `code/08.1-WGBS-TET3-DMRs-whole-genome-computation.R`

## Final tables

- DMR summary/statistics tables and associated processed data are under `data/`.
