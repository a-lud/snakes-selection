# Snake Selection Paper

This repository contains the code, figures and tables presented in the paper:

> Genomic signals of shifts in selection associated with the recent land-sea transition in front-fanged snakes (Elapidae)

The custom pipelines used by some of these scripts can be found at [nf-pipelines][nfpipelines] (build: 2f89cd6) -
specifically the [MSA][msa] and [HyPhy][hyphy] pipelines.

## Directories in the repository

### Manuscript

This directory contains all the documents/figures/tables relating to the published
manuscript.

### Scripts

Contains numbered directories which indicate the order that they were used (for
the most part). 

- **00_manuscript-figures-tables**: Code to generate manuscript materials
- **01_analysis-preparation**: Scripts to build utility data objects
- **02_data-generation-scripts**: Series of directories with scripts used to process the data
- **03_selection**: Code to convert HyPhy output to dataframes and explore selection results
- **04_GO-analysis**: Script to run _topGO_
- **util**: Utility scripts I wrote to help with data processing

[nfpipelines]: https://github.com/a-lud/nf-pipelines
[msa]: https://github.com/a-lud/nf-pipelines/wiki/Multiple-Sequence-Alignment
[hyphy]: https://github.com/a-lud/nf-pipelines/wiki/HyPhy
