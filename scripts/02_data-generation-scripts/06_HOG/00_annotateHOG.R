# ---------------------------------------------------------------------------- #
# Annotate HOGs using OMA annotation
suppressPackageStartupMessages({
  library(tidyverse)
})

# ---------------------------------------------------------------------------- #
# Data
annotation <- read_rds(
  file = 'data/utility-data/master-annotation-16-sample.rds'
)

headers <- read_csv(
  file = 'data/oma-16-samples/hog-headers.csv',
  col_names = TRUE, 
  col_types = cols(),
)

# Gained and lost HOG lists
gained <- list.files(
  path = 'data/HOG-data/hogs-gain-loss',
  pattern = 'blast-gained-samples-threshold-no-hit.txt',
  full.names = TRUE,
  recursive = TRUE
) %>%
  set_names(sub('.+/hogs-gain-loss/(.*)/.*', '\\1', .)) %>%
  map(read_lines)

lost <- list.files(
  path = 'data/HOG-data/hogs-gain-loss',
  pattern = 'blast-lost-samples-no-hit.txt', 
  full.names = TRUE,
  recursive = TRUE
) %>%
  set_names(sub('.+/hogs-gain-loss/(.*)/.*', '\\1', .)) %>%
  map(read_lines)

# Snakes under assessment
samples <- read_lines('data/HOG-data/hogs-gain-loss/ancestral-node-0/analysis-node.txt') %>%
  str_split(string = ., pattern = '/') %>%
  unlist()

# ---------------------------------------------------------------------------- #
# OMA group identifiers
ogids <- read_csv(file = 'data/oma-16-samples/Output/OrthologousGroups.csv') %>%
  pivot_longer(names_to = 'sample', values_to = 'id', cols = 2:length(.)) %>%
  filter(!is.na(id))

# ---------------------------------------------------------------------------- #
# Iterate over conditions - join HOG information with OG inforamtion, then annotations
annotated <- map2(gained[1], lost[1], ~{
  
  # Get only HOGs relating to current analysis
  df_headers_gain <- headers %>%
    filter(HOG %in% .x) %>%
    select(HOG, id = id_clean, sample = snake)
  
  df_headers_lost <- headers %>%
    filter(HOG %in% .y) %>%
    select(HOG, id = id_clean, sample = snake)
  
  # Join the OGIDs with the annotation object - sequence headers + annotation
  df_anno <- ogids %>%
    left_join(annotation)
  
  # Join HOG df with annotation on headers
  df_gain <- left_join(df_headers_gain, df_anno) %>%
    distinct()
  
  df_lost <- left_join(df_headers_lost, df_anno) %>%
    distinct()
  
  return(list('gained' = df_gain, 
              'lost' = df_lost))
})

# ---------------------------------------------------------------------------- #
# Write annotations to file
iwalk(annotated, ~{
  annotated$`ancestral-node-0`$gained
  annotated$`ancestral-node-0`$lost
  
  path_gain <- file.path('data/HOG-data/hogs-gain-loss', .y, 'annotation-gained-hogs.csv')
  path_lost <- file.path('data/HOG-data/hogs-gain-loss', .y, 'annotation-lost-hogs.csv')
  
  write_csv(
    x = .x[['gained']], 
    file = path_gain, 
    col_names = TRUE,
  )
  
  write_csv(
    x = .x[['lost']],
    file = path_lost, 
    col_names = TRUE
  )
})
