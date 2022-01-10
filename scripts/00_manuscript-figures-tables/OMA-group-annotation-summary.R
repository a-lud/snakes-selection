# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
})

# ---------------------------------------------------------------------------- #
# Data
annotation <- read_rds(here('data', 'utility-data', 'master-annotation.rds'))
oma.ids <- fs::dir_ls(
  path = here('data', 'oma-16-samples', 'complete-oma-groups', 'four-sample', 'protein'),
  glob = '*.fa'
) %>%
  sub('.*/(.*).fa', '\\1', .)

# ---------------------------------------------------------------------------- #
# Number of annotated OGIDs
unannotated <- annotation[annotation$Group %in% oma.ids, ] %>%
  ungroup() %>%
  filter(
    Group == Gene,
    across(
      .cols = 3:6,
      .fns = ~is.na(.x)
    )
  ) %>%
  distinct(Group)

# ---------------------------------------------------------------------------- #
# Annotation results
length(oma.ids) - nrow(unannotated) # 3256 annotated genes
100 - ((nrow(unannotated)/length(oma.ids)) * 100) # 90.6% with some annotation
