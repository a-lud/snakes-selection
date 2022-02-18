# ---------------------------------------------------------------------------- #
# Overlap in marine accelerated genes
#
# This script is simply intersecting the 379 genes we identified in A. laevis
# with marine-associated genes reported by the literature.
#
# The process is:
#   1. Get all synonyms for our genes from BiomaRt
#   2. Convert the table to long form (retaining original annotations as key)
#   3. Find matches based on exact gene symbols
#   4. Find matches on remaining synonym matches - ,manually check they're real
#
# This script is to find overlapping genes with the literature that are A. 
# laevis specific!
#
# NOTE: All commented out code has bee run previously, with outputs created
#       to simplify the script.
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
})

fs::dir_create(
  path = here('manuscript', 'figures')
)

# ---------------------------------------------------------------------------- #
# Sea snake data
# marine <- read_lines(
#   here('data', 'selection-results', 'four', 'absrel-pval-0.001-marine.txt')
# )
# 
# marine <- read_csv(
#   file = here('data', 'utility-data', 'master-GO-annotation.csv.gz'),
#   col_names = TRUE,
#   col_type = cols()
# ) %>%
#   filter(
#     Group %in% marine
#   ) %>%
#   pull('Gene') %>%
#   unique()

# ---------------------------------------------------------------------------- #
# Get symbol synonyms for our genes
# ensembl <- biomaRt::useEnsembl(
#   biomart = 'genes',
#   dataset = 'hsapiens_gene_ensembl'
# )
# 
# # Attributes to return
# attributes <- c('external_synonym', "hgnc_symbol")
# 
# # Query to BiomaRt
# marine_synonyms <- biomaRt::getBM(
#   attributes = attributes,
#   filters = "hgnc_symbol",
#   values = marine,
#   mart = ensembl
# ) %>%
#   as_tibble()
# 
# # Add in query sequences as one of the synonyms - synonym column is what we'll use
# # to search
# temp <- tibble(
#   hgnc_symbol = marine,
#   external_synonym = marine
# )
# 
# # Bind biomaRt tibble with all names
# marine_synonyms <- bind_rows(
#   marine_synonyms, 
#   temp
# ) %>%
#   filter(
#     external_synonym != '',
#     ! str_detect(string = hgnc_symbol, pattern = '^OG')
#   ) %>%
#   select(
#     symbol = hgnc_symbol,
#     synonym = external_synonym
#   ) %>%
#   arrange(symbol)
# 
# write_rds(
#   x = marine_synonyms,
#   file = here('data', 'literature-marine-genes', 'marine-psgs.rds')
# )

# ---------------------------------------------------------------------------- #
# Published 

# Chikina et al. 2016
# chikina_690 <- readxl::read_xlsx(
#   path = here('data', 'literature-marine-genes', 'chikina-2016', 'msw112_supplementary_data', 'Supp_Tables_Legends.xlsx'), 
#   sheet = 1, 
#   col_names = TRUE, 
#   skip = 2
# ) %>%
#   filter(
#     `P-value` <= 0.05
#   ) %>%
#   pull(Gene)
# 
# # Foote et al. 2015
# foote_274 <- read_csv(
#   file = here('data', 'literature-marine-genes', 'foote-2015', 'foote_significant.csv'),
#   col_names = TRUE,
#   col_types = cols()
# ) %>%
#   pull(Gene) %>% 
#   unique()
# 
# # Gayk et al. 2018
# Gayk_152 <- readxl::read_xlsx(
#   path = here(
#     'data', 'literature-marine-genes',
#     'gayk-2018', '12862_2018_1181_MOESM3_ESM.xlsx'
#   ), 
#   sheet = 2, 
#   col_names = TRUE
# ) %>% 
#   arrange(-LRT) %>%
#   mutate(
#     Gene_ID = sub('^[^\\.]*\\.', '', Gene_ID)
#   ) %>% 
#   pull(Gene_ID)
# 
# # Li et al. 2021 - Don't actually provide any gene list OR annotation
# # Only PSG I can get are from the GO enrichment results
# # readxl::read_xlsx(
# #   path = here('data', 'literature-marine-genes', 'li-2021', 
# #               'msab212_supplementary_data','revised_Supplementary_Tables.xlsx'),
# #   sheet = 13,
# #   col_names = TRUE
# # ) %>%
# #   select(8) %>%
# #   pull() %>%
# #   str_split(string = ., pattern = ';') %>%
# #   unlist() %>%
# #   sub('\\s+', '', .) %>%
# #   unique()
# 
# # McGowan et al. 2012 - Dolphins
# mcgowan_228 <- readxl::read_xls(
#   path = here('data', 'literature-marine-genes', 'mcgowan-2012', 'rspb20120869supp3.xls'),
#   sheet = 1,
#   col_names = TRUE
# ) %>%
#   filter(`dN/dS` > 1) %>%
#   select(`Associated Gene Name`) %>%
#   pull()
# 
# # Peng et al. 2020
# peng_470 <- read_csv(
#   file = here(
#     'data', 'literature-marine-genes', 
#     'peng-2020', 'tabula-rough-export-from-pdf.csv'
#   ), 
#   col_names = c('custom_id', 'symbol', 'description'),
#   col_types = cols(), 
#   na = ""
# ) %>%
#   select(symbol) %>%
#   filter(
#     !is.na(symbol),
#     nchar(symbol) > 2,
#     ! symbol == '0278.1|'
#   ) %>%
#   mutate(
#     symbol = str_replace(
#       string = symbol, 
#       pattern = 'gi\\|602626649\\|ref\\|XP_00742', 
#       replacement = 'gi\\|602626649\\|ref\\|XP_007420278.1'
#     ),
#     symbol = str_replace(
#       string = symbol,
#       pattern = 'ETE71823\\.',
#       replacement = 'ETE71823.1'
#     ),
#     symbol = sub('.*\\|', '', symbol),
#     symbol = sub('_.*', '', symbol)
#   ) %>%
#   filter(symbol != "XP") %>%
#   pull(symbol)
# 
# # Sun et al. 2012 - Dolphin Ensembl IDs do not match Ensembl V66 - can't
# # actually convert to a useful identifier...
# # sun_693 <- readxl::read_xls(
# #   path = here(
# #     'data', 'literature-marine-genes','sun-2012',
# #     'evs123_Supplementary_Data',
# #     'Supplementary_Table_2_PSGs_and_Functinal_clusters.xls'
# #   ),
# #   sheet = 1,
# #   col_names = TRUE
# # ) %>%
# #   filter(FDR <= 0.05) %>%
# #   pull(`Ensembl GeneID`)
# #
# # dolphin_gtf <- rtracklayer::import(
# #   '~/Desktop/Tursiops_truncatus.turTru1.66.gtf.gz', 
# # ) %>%
# #   as.data.frame() %>%
# #   as_tibble() %>%
# #   select(gene_id, gene_name)

# Store each dataset as list
# lit_genes <- list(
#   'chikina' = chikina_690,
#   'foote' = foote_274,
#   'gayk' = Gayk_152,
#   'mcgowan' = mcgowan_228,
#   'peng' = peng_470
# )

# write_rds(
#   x = lit_genes,
#   file = 'data/literature-marine-genes/list-literature-genes.rds',
#   compress = 'gz'
# )

# ---------------------------------------------------------------------------- #
# Get gene synonyms from biomaRt

# Rerun ONLY if need to regenerate the output file - takes a while otherwise
# gene_synonyms <- imap_dfr(lit_genes, ~{
#   df <- biomaRt::getBM(
#     attributes = attributes,
#     filters = "hgnc_symbol",
#     values = .x,
#     mart = ensembl
#   ) %>%
#     as_tibble()
# 
#   # Append original query symbols to synonym column - this column is used for searching
#   df2 <-
#     tibble(
#       hgnc_symbol = .x,
#       external_synonym = .x
#     )
# 
#   # Join query tibble with original symbol tibble - filter out junk
#   df_bind <- bind_rows(
#     df, df2
#   ) %>%
#     filter(external_synonym != '') %>%
#     select(
#       symbol = hgnc_symbol,
#       synonym = external_synonym
#     ) %>%
#     arrange(symbol) %>%
#     mutate(
#       dataset = .y
#     )
# }) %>%
#   distinct() %>%
#   mutate(
#     symbol = toupper(symbol),
#     synonym = toupper(synonym),
#     dataset = case_when(
#       dataset == 'chikina' ~ 'Chikina et al. 2016\n(Marine mammals)',
#       dataset == 'foote' ~ 'Foote et al. 2015\n(Marine mammals)',
#       dataset == 'gayk' ~ 'Gayk et al. 2018\n(Diving bird)',
#       dataset == 'mcgowan' ~ 'McGowan et al. 2012\n(Dolphin)',
#       dataset == 'peng' ~ 'Peng et al. 2020\n(Sea snake)'
#     )
#   )

# write_csv(
#   x = gene_synonyms,
#   file = here('data', 'literature-marine-genes', 'literature-marine-genes.csv')
# )

# ---------------------------------------------------------------------------- #
# Join tables to find overlap
# marine_synonyms <- read_rds(here('data', 'literature-marine-genes', 'marine-psgs.rds'))
# 
# gene_synonyms <- read_csv(
#   file = here('data', 'literature-marine-genes', 'literature-marine-genes.csv'),
#   col_names = TRUE,
#   col_types = cols()
# )
# 
# # 1. Get matches where the annotations are the same
# match_17 <- unique(marine_synonyms$symbol)[unique(marine_synonyms$symbol) %in% unique(gene_synonyms$symbol)]
# 
# # 2. Get matches based on synonym - remove the matches above (as they're accounted for) - and write
# #    the file to be checked manually
# marine_synonym <- marine_synonyms %>% filter(! symbol %in% match_17)
# lit_synonym <- gene_synonyms %>% filter( ! gene_synonyms$synonym %in% match_17 )
# 
# overlap_synonym <- inner_join(marine_synonym, lit_synonym, by = 'synonym') %>% arrange(symbol.x)
# 
# write_csv(
#   x = overlap_synonym,
#   file = here(
#     'data', 'literature-marine-genes', 'synonym-overlap-379-check.csv'
#   ),
#   col_names = TRUE
# )
# 
# # After manual curation, none of the synonymous matches are legitimate. The shared synonyms
# # are due to different genes having the same previous name (symbol) - they're now different
# overlap_379 <- gene_synonyms %>% filter(symbol %in% match_17) %>%
#   select(symbol, dataset) %>%
#   distinct() %T>%
#   write_csv(
#     file = here(
#       'manuscript', 'tables',
#       'supplementary-psg-overlap-literature.csv'
#     ), 
#     col_names = TRUE
#   )
# 
# write_csv(
#   x = overlap_379,
#   file = here('data/literature-marine-genes/overlap-psgs-literature.csv'),
#   col_names = TRUE
# )

overlap_379 <- read_csv(
  file = here('data', 'literature-marine-genes', 'overlap-psgs-literature.csv'),
  col_names = TRUE,
  col_types = cols()
) %>%
  group_by(dataset) %>%
  summarise(count = n()) %>%
  mutate(
    n_genes = case_when(
      dataset == 'Chikina et al. 2016\n(Marine mammals)' ~ '690',
      dataset == 'Foote et al. 2015\n(Marine mammals)' ~ '274',
      dataset == 'Gayk et al. 2018\n(Diving bird)' ~ '152',
      dataset == 'McGowan et al. 2012\n(Dolphin)' ~ '228',
      dataset == 'Peng et al. 2020\n(Sea snake)'~ '470'
    ),
    dataset = factor(x = dataset, levels = c(unique(dataset), 'McGowan et al. 2012\n(Dolphin)'))
  )

# ragg::agg_png(
#   filename = here('manuscript', 'figures', 'supplementary-literature-overlap.png'),
#   width = 1500,
#   height = 1500,
#   units = 'px',
#   res = 144
# )
p <- overlap_379 %>%
  ggplot(
    aes(
      x = dataset,
      y = count,
      fill = dataset
    )
  ) +
  geom_bar(stat = 'identity', colour = 'black', alpha = 0.8) +
  geom_text(aes(label = n_genes), vjust = -0.5, size = 5) +
  ggpomological::scale_fill_pomological() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = seq(1,10, 1), limits = c(0, 10)) +
  labs(
    y = 'Number of Genes'
  ) +
  theme_bw() +
  theme(
    # Axis labels
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.title.x = element_blank(),
    
    # Axis text
    axis.text = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Legend
    # legend.title = element_blank(),
    legend.position = 'none',
    # legend.text = element_text(size = 14)
  )

write_rds(
  x = p, 
  file = here('data/literature-marine-genes/supplementary-psg-lit-overlap.rds')
)
# invisible(dev.off())

