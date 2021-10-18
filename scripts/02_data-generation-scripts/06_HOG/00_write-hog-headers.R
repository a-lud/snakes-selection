# --------------------------------------------------------------------------- #
# Write headers to file for each snake
#
# The headers files created by this script are to be used by the RNA-alignments
# script to subset the bams for HOGs of interest.
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
})

# Samples that can use cleaned identifiers
clean_vec <- c('aipysurusLaevis_ACAGTG_ACTTGA_CGATGT', 'aipysurusLaevis',
               'aipysurusLaevis_KLS0468', 'aipysurusLaevis_L1_eye')

# Samples that must use original
original_vec <- c('aipysurusMosaicus', 'Atenuis_CAGATC_TTAGGC_GCCAAT',
                  'hydrelapsDarwiniensis', 'hydrophisKingii', 'hydrophisMajor',
                  'hydrophisMajor_KLS0460', 'hydrophisPeronii')

# --------------------------------------------------------------------------- #
# Data import

# HOGs mapped to Headers
headers <- vroom::vroom(
  file = 'data/oma-16-samples/hog-headers.csv',
  delim = ',',
  col_names = TRUE,
  col_types = cols()
)

# HOG identifiers
hog_ids <- list.files(
  path = 'data/HOG-data/hogs-gain-loss',
  pattern = 'blast.*.txt',
  full.names = TRUE,
  recursive = TRUE
) %>%
  grep(pattern = 'blast-gained-samples-threshold-no-hit.txt', x = ., value = TRUE) %>%
  set_names(sub(
    '.+/hogs-gain-loss/(.*)/blast-gained-samples-threshold-no-hit.txt', 
    '\\1', 
    .)) %>%
  map(read_lines)

# --------------------------------------------------------------------------- #
# Write to file by dataset
iwalk(hog_ids, ~{
  df_split <- headers %>%
    filter(HOG %in% .x) %>%
    select(HOG, snake, id, id_clean) %>%
    distinct() %>%
    split(x = ., .$snake)
  
  # Write identifiers for each sample to file
  walk(names(df_split), function(snake) {
    path <- file.path(
      'data/HOG-data/hogs-gain-loss',
      .y,
      'hog-headers-by-snake'
    )
    dir.create(path = path, showWarnings = FALSE, recursive = TRUE)
    
    file <- file.path(
      path,
      paste0(snake, '.txt')
    )
    
    if(snake %in% clean_vec) {
      write_lines(
        x = df_split[[snake]][['id']],
        file = file
      )
    } else if(snake %in% original_vec) {
      # Multiply the length by 3 to get CDS length
      ids <- df_split[[snake]] %>%
        mutate(
          num = sub('.*:1-(.*)', '\\1', id),
          num = as.integer(num),
          num = num * 3,
          id_cds = paste0(id_clean, ":1-", num)
        ) %>%
        extract2('id_cds')
      
      write_lines(
        x = ids,
        file = file
      )
    }
  })
})
