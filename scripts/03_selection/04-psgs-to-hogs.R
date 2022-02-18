# ---------------------------------------------------------------------------- #
# Get HOGs relating to PSGs

# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(furrr)
  library(here)
})

future::plan("multisession", workers = 4)
fs::dir_create(path = here('data', 'psg-hogs'))

# ---------------------------------------------------------------------------- #
# OMA-groups for PSGs
# smpl <- c('Group', 'aipysurusLaevis', 
#              'notechisScutatus', 'pseudonajaTextilis', 'najaNaja')

# OMA-group data
all.oma <- read_csv(
  file = here('data', 'oma-16-samples', 'Output', 'OrthologousGroups.csv'),
  col_names = TRUE,
  col_types = cols()
) %>%
  select(c('Group', 'aipysurusLaevis'))

psgs <- read_lines(
  here('data', 'selection-results', 'four', 'absrel-pval-0.001-marine.txt')
)

psgs.oma <- all.oma %>% 
  filter(Group %in% psgs) %>%
  pivot_longer(names_to = 'samples', values_to = 'ID', 2)

# ---------------------------------------------------------------------------- #
# HOG data
hogs <- fs::dir_ls(
  path = here('data', 'oma-16-samples', 'Output', 'HOGFasta'),
  glob = "*.fa", 
) %>%
  set_names(value = sub('.*/(.*).fa', '\\1', .)) %>%
  future_imap_dfr(.progress = TRUE, ~{
    lines <- read_lines(.x) %>%
      grep(pattern = "^>", x = ., value = TRUE)
    tibble(
      "HOG" = .y,
      "ID" = lines
    )
  }) %>%
  separate(
    col = ID, into = c('ID', 'samples'),
    sep = ' '
  ) %>%
  mutate(
    ID = sub('^>', '', ID),
    samples = gsub('\\[|\\]', '', samples)
  )

# ---------------------------------------------------------------------------- #
# Join OMA-group with HOG data using sequence identifiers + species as keys
psgs.omag.hog <- left_join(psgs.oma, hogs)

psg.hogs <- psgs.omag.hog %>% pull(HOG) %>% unique()


# ---------------------------------------------------------------------------- #
# Write HOGs to file
write_lines(
  x = psg.hogs, 
  file = here('data', 'psg-hogs', 'psg-hog-list.txt')
)
