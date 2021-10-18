# ---------------------------------------------------------------------------- #
# Summaries to GO
# 
# Taken from Steve Pederson's work from:
# https://uofabioinformaticshub.github.io/summaries2GO/MakeSummaries

# Libraries
suppressPackageStartupMessages({
  library(GO.db)
  library(graph)
  library(dnet)
  library(magrittr)
  library(tidyverse)
})

# ---------------------------------------------------------------------------- #
# GO Summaries - pruning tree

# Loading Steve's pre-made
goSummaries <- url("https://uofabioinformaticshub.github.io/summaries2GO/data/goSummaries.RDS") %>%
  readRDS()

# Making up-to-date goSummaries manually (Steve's code)
graphs <- c(BP = "bp", CC = "cc", MF = "mf") %>%
  lapply(makeGOGraph) %>%
  lapply(function(x){removeNode("all", x)}) %>%
  lapply(dDAGreverse)
write_rds(x = graphs, file = 'data/go-graph.rds', compress = 'gz')

currentGoSummaries <- lapply(graphs, function(x){
  lng <- dDAGlevel(x, "longest_path") - 1
  shrt <- dDAGlevel(x, "shortest_path") - 1
  tips <- dDAGtip(x)
  tibble(
    id = unique(c(names(lng), names(shrt))),
    shortest_path = shrt,
    longest_path = lng,
    terminal_node = id %in% tips
  )
}) %>%
  bind_rows() %>%
  mutate(ontology = Ontology(id))
write_rds(x = currentGoSummaries, file = 'data/updatedGOSummary.rds', compress = 'gz')

# ---------------------------------------------------------------------------- #
# Gene length information
lst_geneLength <- list.files(path = 'data/msa', 
                             pattern = 'stats.txt',
                             full.names = TRUE,
                             recursive = TRUE) %>%
  set_names(sub('.+/(.*)/stats.txt', '\\1', .)) %>%
  map(read_tsv, col_names = TRUE, col_types = cols()) %>%
  map(~{
    .x %>%
      select(ogid = file, length = min_len) %>%
      mutate(ogid = sub('.+/(.*)_translated.fasta', '\\1', ogid),
             length = as.integer(length))
  })

walk(names(lst_geneLength), function(analysis){
  write_csv(x = lst_geneLength[[analysis]], 
            file = file.path(paste0('data/msa/', analysis, '/'), 
                             paste0(analysis, '-gene-lengths.csv')), 
            col_names = TRUE)
})
