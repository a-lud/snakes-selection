# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
})

# ---------------------------------------------------------------------------- #
# Import HyPhy data
bad <- c('VRK1', 'GLA') # These were quite significant in N. naja so probs leave

hyphy_d4 <- read_rds(
  file = here('data', 'hyphy-dataframes', 'hyphy-four-sample.rds')
  # file = 'data/hyphy-dataframes/hyphy-four-sample.rds'
)

hyphy_d9 <- read_rds(
  file = here('data', 'hyphy-dataframes', 'hyphy-nine-sample.rds')
)

hyphy <- list(
  'Four' = hyphy_d4,
  'Nine' = hyphy_d9
)

rm(hyphy_d4)
rm(hyphy_d9)

gc()

genes_marine_D4 <- read_lines(here('data', 'selection-results', 'four', 'absrel-pval-0.001-marine.txt'))
genes_marine_D9 <- read_lines(here('data', 'selection-results', 'nine', 'absrel-pval-0.001-marine.txt'))

genes_ter_D4 <- read_lines(here('data', 'selection-results', 'four', 'absrel-pval-0.001-terrestrial.txt'))
genes_ter_D9 <- read_lines(here('data', 'selection-results', 'nine', 'absrel-pval-0.001-terrestrial.txt'))

annotation <- read_csv(
  file = here('data', 'utility-data', 'master-GO-annotation.csv'),
  col_names = TRUE, 
  col_types = cols()
)

gosummaries <- read_rds(
  file = here('data', 'utility-data', 'updatedGOSummary.rds')
) %>%
  filter(
    terminal_node == TRUE,
    shortest_path >= 3,
    longest_path >= 3
  )

# Create summary tables
absrel_sig <- imap_dfr(hyphy, .id = 'Dataset', ~{
  if (.y == 'Four') {
    species <- 'aipysurusLaevis'
    genes_marine <- genes_marine_D4
    # genes_ter <- genes_ter_D4
  } else {
    species <- 'Node_aipysurusLaevis'
    genes_marine <- genes_marine_D9
    # genes_ter <- genes_ter_D9
  }
  
  df <- .x[['absrel']][['Branch Attributes']] %>%
    filter(
      `Corrected P-value` <= 0.001,
      Species == species
    ) %>%
    separate(col = ID, into = c('Model', 'OGID'), sep = '-') %>%
    select(
      OGID, Species, 
      `P-value (aBSREL)` = `Corrected P-value`
    ) %>%
    filter(
      OGID %in% genes_marine
    )
})

relax_sig <- imap_dfr(hyphy, .id = 'Dataset', ~{
  if (.y == 'Four') {
    species <- 'aipysurusLaevis'
    genes_marine <- genes_marine_D4
    # genes_ter <- genes_ter_D4
  } else {
    species <- 'Node_aipysurusLaevis'
    genes_marine <- genes_marine_D9
    # genes_ter <- genes_ter_D9
  }
  
  df <- .x[['relax']][['Test Results']] %>%
    filter(`p-value` <= 0.001) %>%
    rename(
      `K-value` = `relaxation or intensification parameter`,
      `P-value (RELAX)` = `p-value`
    ) %>%
    separate(col = ID, into = c('Model', 'OGID'), sep = '-') %>%
    mutate(
      `Selection` = ifelse(`K-value` > 1,
                           'Intensification',
                           'Relaxation')
    ) %>%
    select(OGID, `P-value (RELAX)`, `K-value`, Selection) %>%
    filter(OGID %in% genes_marine)
})

# Joining objects
absrel_relax <- full_join(
  absrel_sig, relax_sig
) %>%
  left_join(
    annotation,
    by = c('OGID' = 'Group')
  ) %>%
  distinct() %>%
  select(
    Dataset, Gene, Species, `P-value (aBSREL)`,
    `P-value (RELAX)`, `K-value`,
    Selection, `GO-term` = GO, Term = TERM,
    Definition = DEFINITION, Ontology = ONTOLOGY
  ) %>%
  filter(
    ! Gene %in% bad
  )

# Eh - think we drop this
# meme_sig <- hyphy$Nine$meme$MLE %>%
#   separate(col = ID, into = c('Model', 'OGID'), sep = '-') %>%
#   group_by(OGID) %>%
#   mutate(
#     Position = 1:n()
#   ) %>%
#   group_map(.keep = TRUE , .f = ~{
#     .x %>%
#       arrange(`p-value`) %>%
#       mutate(
#         `Adj. P-value` = p.adjust(p = `p-value`, method = 'bonferroni')
#       )
#   }) %>%
#   bind_rows() %>%
#   filter(`p-value` <= 0.1) %>%
#   left_join(
#     annotation,
#     by = c('OGID' = 'Group')
#   ) %>%
#   filter(
#     OGID %in% genes_marine_D9,
#     ! Gene %in% c('GLA', 'VRK1')
#   ) %>%
#   select(
#     Gene, Position, `p-value`, `Adj. P-value`
#   ) %>%
#   distinct()


# Building pretty tables - aBSREL/RELAX table - Annotation table
# Facet by dataset
fs::dir_create(
  path = here('data', 'selection-results', 'summary-tables'),
  recurse = TRUE
)

# General summary
general_summary <- absrel_relax %>%
  select(1:7) %>%
  distinct() %>%
  split(x = ., f = .$Dataset, drop = TRUE) %T>%
  iwalk(~{
    path <- file.path(
      here(
        'data',
        'selection-results',
        'summary-tables',
        paste0(.y, '-summary.csv')
      )
    )
    
    write_csv(
      x = .x,
      file = path,
      col_names = TRUE
    )
  })

# GO term summary
general_go_term_summary <- absrel_relax %>%
  select(1:3, 8:11) %>%
  filter(
    across(
      .cols = 4:7,
      .fns = ~!is.na(.x)
    )
  ) %>%
  distinct() %>%
  split(., .$Dataset, drop = TRUE) %T>%
  iwalk(~{
    path <- file.path(
      here(
        'data',
        'selection-results',
        'summary-tables',
        paste0(.y, '-annotation.csv') 
      )
    )
    
    write_csv(
      x = .x,
      file = path,
      col_names = TRUE
    )
  })

# filtered annotation
filtered_go_term_summary <- absrel_relax_anno %>%
  map(~{
    .x %>%
      filter(`GO-term` %in% gosummaries$id)
  }) %T>%
  imap(~{
    path <- file.path(
      here(
        'data',
        'selection-results',
        'summary-tables',
        paste0(.y, '-annotation-filtered.csv') 
      )
    )
    .x %>%
      write_csv(
        file = path,
        col_names = TRUE
      )
  })

# Write single XLSX file with all tables in it
l <- list(
  'Selection-tests-D4' = general_summary$Four,
  'GO-Summary-D4' = general_go_term_summary$Four,
  'Filtered-GO-summary-D4' = filtered_go_term_summary$Four,
  'Selection-tests-D9' = general_summary$Nine,
  'GO-Summary-D9' = general_go_term_summary$Nine,
  'Filtered-GO-summary-D9' = filtered_go_term_summary$Nine
)

openxlsx::write.xlsx(
  x = l, 
  file = here(
    'data',
    'selection-results',
    'summary-tables',
    'summary-table.xlsx'
  )
)


