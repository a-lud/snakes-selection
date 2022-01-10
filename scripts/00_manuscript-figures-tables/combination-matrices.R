# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
})

# ---------------------------------------------------------------------------- #
# Load combination matrix rds object
comb <- read_rds(
  file = here('data', 'selection-results', 'combination-matrices.rds')
  )

comb_four <- comb$four$`001`
comb_nine <- comb$nine$`001`

# ---------------------------------------------------------------------------- #
# Write combination matrix to file
tibble(
  id = names(ComplexHeatmap::comb_size(comb_four)),
  intersect = ComplexHeatmap::comb_size(comb_four)
) %>%
  separate(
    col = id,
    into = c('Aipysurus laevis', 'Notechis scutatus', "Pseudonaja textilis", "Naja naja"),
    sep = "(?<=.)"
  ) %>%
  arrange(-intersect) %>%
  write_csv(
    file = here('manuscript', 'tables', 'combination-matrix.csv'),
    col_names = TRUE
  )

# tibble(
#   id = names(ComplexHeatmap::comb_size(comb_nine)),
#   intersect = ComplexHeatmap::comb_size(comb_nine)
# ) %>%
#   separate(
#     col = id,
#     into = c('Node Aipysurus laevis', 'Node Notechis scutatus', "Node Pseudonaja textilis", "Naja naja"),
#     sep = "(?<=.)"
#   ) %>%
#   arrange(-intersect) %>%
#   write_csv(
#     file = here('data', 'general-summary-tables', 'combination-mat-nine.csv'),
#     col_names = TRUE
#   )
# 
