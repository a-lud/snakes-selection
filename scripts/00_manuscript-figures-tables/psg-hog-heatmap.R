# ---------------------------------------------------------------------------- #
# Heatmap: n-seqs in each HOG by snake

# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
})

# ---------------------------------------------------------------------------- #
# Data
smpl <- c('aipysurusLaevis', 'notechisScutatus', 'pseudonajaTextilis', 'najaNaja')

hogs <- read_lines(
  here('data', 'psg-hogs', 'psg-hog-list.txt')
)

phyletic.profile <- read_tsv(
  file = here('data', 'oma-16-samples', 'Output', 'PhyleticProfileHOGs.txt'),
  col_names = TRUE, 
  col_types = cols(),
  skip = 4
)

# ---------------------------------------------------------------------------- #
# Subset for HOGs of interest
phyletic.profile %>%
  filter(Group %in% hogs) %>%
  select(all_of(c('Group', smpl))) %>%
  # column_to_rownames('Group') %>%
  # pheatmap::pheatmap()
  pivot_longer(names_to = 'samples', values_to = 'ortholog_counts', 2:5) %>%
  tidyHeatmap::heatmap(
    .row = Group, 
    .column = samples, 
    .value = ortholog_counts, 
    .scale = "none",
    palette_value = viridisLite::viridis(n = 3)
  )
