# ---------------------------------------------------------------------------- #
# Transcriptome cleanup: Subset GFF files for HiQ sequences only
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Libraries
library(readr)
library(tibble)
library(purrr)
library(magrittr)

# ---------------------------------------------------------------------------- #
# Data
hiq_headers <- list.files(
  path = '../../analyses/genome-data/transcriptomes/completeORF',
  pattern = '.cds',
  full.names = TRUE
) %>%
  set_names(sub('.completeORF.cds', '', basename(.))) %>%
  map(~{
    lines <- read_lines(file = .x) %>%
      grep(pattern = '>', x = ., value = TRUE) %>%
      sub('>(.*)\\..*', '\\1', .)
  })

gffs <- list.files(
  path = '../../analyses/genome-data/transcriptomes/transdecoder',
  pattern = '.gff3',
  full.names = TRUE
) %>%
  set_names(sub('.cdhit.fasta.transdecoder.gff3', '', basename(.))) %>%
  map(~{
    rtracklayer::readGFF(filepath = .x, version = 3) %>%
      as.data.frame() %>%
      as_tibble()
  })

# ---------------------------------------------------------------------------- #
# Subset GFF
gffs_filtered <- map(names(gffs), function(sample){
  h <- hiq_headers[[sample]]
  g <- gffs[[sample]]
  
  g <- dplyr::filter(g, seqid %in% h)
}) %>%
  set_names(names(gffs))

# ---------------------------------------------------------------------------- #
# Free up some space
rm(gffs)
gc()

# ---------------------------------------------------------------------------- #
# Write GFFs to file
walk(names(gffs_filtered), function(name){
  rtracklayer::export.gff(object = gffs_filtered[[name]], paste0('../../analyses/genome-data/transcriptomes/completeORF/', name, '.gff3'))
})
