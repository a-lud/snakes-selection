# Snake Selection Paper

This repository contains scripts used for selection and HOG analyses. In addition to
these scripts, the repository [nf-pipelines][nfpipelines]) was also used to generate data -
specifically the [MSA][msa]) and [HyPhy][hyphy]) pipelines.

## Structure of the repository

There are three main directories within this repository that house most of the information relating
to this project:

- **data**: Contains most of the results in more user-friendly formats. Specifically, files relating
  to HyPhy results, selection-summary results along with utility data (annotaiton object, GO annotations)
- **scripts**: This directory contains all the scripts used to generate the auxiliary data files,
  run the analyses and explore the results. These should all be numbered in the order they were used.
- **figures**: Figures generated form the exploration of the data. Should house all figures used
  in the publication/some extra ones.

[nfpipelines]: https://github.com/a-lud/nf-pipelines
[msa]: https://github.com/a-lud/nf-pipelines/wiki/Multiple-Sequence-Alignment
[hyphy]: https://github.com/a-lud/nf-pipelines/wiki/HyPhy
