suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
  library(patchwork)
})

# ---------------------------------------------------------------------------- #
# A. laevis PSGs that intersect with literature genes
psg.overlap.lit.17 <- read_rds(here('data/literature-marine-genes/supplementary-psg-lit-overlap.rds'))
non.sig.overlap.lit.154 <-  read_rds(here('data/literature-marine-genes/supplementary-non-psg-lit-overlap.rds'))


p <- psg.overlap.lit.17/non.sig.overlap.lit.154 + plot_annotation(tag_levels = 'A')

ragg::agg_png(
  filename = here('manuscript', 'figures', 'literature-overlap-faceted.png'),
  width = 1500,
  height = 1500,
  units = 'px',
  res = 144
)
print(p)
invisible(dev.off())
