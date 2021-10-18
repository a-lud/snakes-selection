#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(magrittr)
})

# --------------------------------------------------------------------------- #
# Arguments
parser <- ArgumentParser(
  description = 'Get HOG sequences by group'
)

parser$add_argument('-f', '--file', help = 'PhyleticProfile path')
parser$add_argument('-s', '--sample', help = 'Sample names')
parser$add_argument('-b', '--background', help = 'Background samples')
parser$add_argument('-o', '--outdir', help = 'Output directory')

args <- parser$parse_args()

# --------------------------------------------------------------------------- #
# Pipeline

# Import phyleticProfile
df <- read_tsv(
  file = args$file,
  col_names = TRUE,
  col_types = cols(),
  skip = 4
)

# Convert group-strings to vectors
args$sample <- unlist(str_split(string = args$sample, pattern = ' '))
args$background <- unlist(str_split(string = args$background, pattern = ' '))

# Wiggle room only applies when N-samples >= 5
nthreshold <- ceiling(length(args$sample)/2)

# Get gained within sample group
unique_gain_sample <- df %>%
  filter(
    if_all(args$sample, ~ .x != 0),
    if_all(args$background, ~ .x == 0)
  ) %>%
  mutate(
    Group = sub('HOG', '', Group),
    Group = sub('^0+', '', Group),
    Group = paste0('HOG', Group)
  ) %>%
  extract2('Group')

# Get lost within sample group
unique_lost_sample <- df %>%
  filter(
    if_all(args$sample, ~ .x == 0),
    if_all(args$background, ~ .x != 0)
  ) %>%
  mutate(
    Group = sub('HOG', '', Group),
    Group = sub('^0+', '', Group),
    Group = paste0('HOG', Group)
  ) %>%
  extract2('Group')

# --------------------------------------------------------------------------- #
# Get HOGs shared by N-samples within group
unique_threshold_gain_sample <- df %>%
  filter(
    if_all(.cols = args$background, .fns = ~.x == 0),
    if_any(.cols = args$sample, .fns = ~.x != 0),
    rowSums(select(., args$sample) > 0) >= nthreshold
  ) %>%
  mutate(
    Group = sub('HOG', '', Group),
    Group = sub('^0+', '', Group),
    Group = paste0('HOG', Group)
  ) %>%
  extract2('Group')


# --------------------------------------------------------------------------- #
# Write to file
write_lines(
  x = unique_gain_sample,
  file = paste0(args$outdir, '/gained-samples-id.txt')
)

write_lines(
  x = unique_lost_sample,
  file = paste0(args$outdir, '/lost-samples-id.txt')
)

write_lines(
  x = unique_threshold_gain_sample,
  file = paste0(args$outdir, '/gained-samples-threshold-id.txt')
)

