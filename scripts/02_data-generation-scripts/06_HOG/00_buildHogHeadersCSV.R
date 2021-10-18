# ---------------------------------------------------------------------------- #
# Build CSV containing all HOG IDs and their respective headers

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(furrr)
})

# For parallel work
plan(multisession, workers = 6)

# --------------------------------------------------------------------------- #
# Read in all HOG fasta files as BSstring objects
files <- list.files(path = 'data/oma-16-samples/Output/HOGFasta', pattern = '.fa', full.names = TRUE) %>%
  set_names(sub('.fa', '', basename(.)))

bs_fasta <- future_map(files, ~{
  Biostrings::readAAStringSet(
    filepath = .x, 
    use.names = TRUE
  )
})

# --------------------------------------------------------------------------- #
# Create clean headers
hog_headers <- future_imap(bs_fasta, ~{
  df <- as_tibble(as.data.frame(.x@ranges))[,3:4] %>%
    separate(col = names, into = c('id', 'snake'), sep = ' ') %>%
    mutate(
      id_clean = sub(':.*', '', id),
      snake = gsub('\\[|\\]', '', snake),
      HOG = .y
    ) %>%
    select(HOG, snake, id_clean, id, width)
}) %>%
  bind_rows()

# --------------------------------------------------------------------------- #
# Output
write_csv(
  x = hog_headers,
  file = 'data/oma-16-samples/hog-headers.csv',
  col_names = TRUE
)
