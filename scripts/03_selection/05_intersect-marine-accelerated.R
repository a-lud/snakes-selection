# ---------------------------------------------------------------------------- #
# Overlap in marine accelerated genes
#
# This script is simply intersecting the 379 genes we identified in A. laevis
# with the 690 significant marine-accelerated genes reported by Chikina et al.
# 2016.
#
# The process is:
#   1. Get all synonyms for our genes from BiomaRt
#   2. Convert the table to long form (retaining original annotations as key)
#   3. Intersect the table from Chikina et al. with ours
#   4. See what overlap there is
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
})

# ---------------------------------------------------------------------------- #
# Sea snake data
marine_379 <- read_lines(
  here('data', 'selection-results', 'four', 'absrel-pval-0.001-marine.txt')
)

marine_379 <- read_csv(
  file = here('data', 'utility-data', 'master-GO-annotation.csv.gz'),
  col_names = TRUE,
  col_type = cols()
) %>%
  filter(
    Group %in% marine_379
  ) %>%
  extract2('Gene') %>%
  unique()

marine_6 <- read_lines(
  here('data', 'selection-results', 'nine', 'absrel-pval-0.001-marine.txt')
)

marine_6 <- read_csv(
  file = here('data', 'utility-data', 'master-GO-annotation.csv.gz'),
  col_names = TRUE,
  col_type = cols()
) %>%
  filter(
    Group %in% marine_6
  ) %>%
  extract2('Gene') %>%
  unique()

# remove two borderline genes
marine_6 <- marine_6[-c(5,6)]

# ---------------------------------------------------------------------------- #
# Get symbol synonyms for our genes
ensembl <- biomaRt::useEnsembl(
  biomart = 'genes',
  dataset = 'hsapiens_gene_ensembl'
)

# Attributes to return
attributes <- c('external_synonym', "hgnc_symbol")

# Query to BiomaRt
marine_379_synonyms <- biomaRt::getBM(
  attributes = attributes,
  filters = "hgnc_symbol",
  values = marine_379,
  mart = ensembl
) %>%
  as_tibble()

# Add in query sequences as one of the synonyms - synonym column is what we'll use
# to search
temp <- tibble(
  hgnc_symbol = marine_379,
  external_synonym = marine_379
)

# Bind biomaRt tibble with all names
marine_379_synonyms <- bind_rows(
  marine_379_synonyms, 
  temp
) %>%
  filter(
    external_synonym != '',
    ! str_detect(string = hgnc_symbol, pattern = '^OG')
  ) %>%
  select(
    symbol = hgnc_symbol,
    synonym = external_synonym
  ) %>%
  arrange(symbol)

# Six genes
marine_6_synonyms <- biomaRt::getBM(
  attributes = attributes,
  filters = "hgnc_symbol",
  values = marine_6,
  mart = ensembl
) %>%
  as_tibble()

temp <- tibble(
  hgnc_symbol = marine_6,
  external_synonym = marine_6
)

marine_6_synonyms <- bind_rows(
  marine_6_synonyms, 
  temp
) %>%
  filter(
    external_synonym != '',
    ! str_detect(string = hgnc_symbol, pattern = '^OG')
  ) %>%
  select(
    symbol = hgnc_symbol,
    synonym = external_synonym
  ) %>%
  arrange(symbol)

# ---------------------------------------------------------------------------- #
# Published 

# Chikina et al. 2016
chikina_690 <- readxl::read_xlsx(
  path = here('data', 'literature-marine-genes', 'chikina-2016', 'msw112_supplementary_data', 'Supp_Tables_Legends.xlsx'), 
  sheet = 1, 
  col_names = TRUE, 
  skip = 2
) %>%
  filter(
    `P-value` <= 0.05
  ) %>%
  pull(Gene)

# Foote et al. 2015
foote_191 <- read_csv(
  file = here('data', 'literature-marine-genes', 'foote-2015', 'edited.csv'),
  col_names = c('locus', 'symbol', 'pval'),
  col_types = cols(), 
  col_select = c('symbol', 'pval')
) %>%
  pull(symbol)

# Gayk et al. 2018
Gayk_152 <- readxl::read_xlsx(
  path = here(
    'data', 'literature-marine-genes',
    'gayk-2018', '12862_2018_1181_MOESM3_ESM.xlsx'
  ), 
  sheet = 2, 
  col_names = TRUE
) %>% 
  arrange(-LRT) %>%
  mutate(
    Gene_ID = sub('^[^\\.]*\\.', '', Gene_ID)
  ) %>% 
  pull(Gene_ID)

# Li et al. 2021 - Don't actually provide any gene list OR annotation
# Only PSG I can get are from the GO enrichment results
# readxl::read_xlsx(
#   path = here('data', 'literature-marine-genes', 'li-2021', 
#               'msab212_supplementary_data','revised_Supplementary_Tables.xlsx'),
#   sheet = 13,
#   col_names = TRUE
# ) %>%
#   select(8) %>%
#   pull() %>%
#   str_split(string = ., pattern = ';') %>%
#   unlist() %>%
#   sub('\\s+', '', .) %>%
#   unique()

# McGowan et al. 2012 - Dolphins
mcgowan_228 <- readxl::read_xls(
  path = here('data', 'literature-marine-genes', 'mcgowan-2012', 'rspb20120869supp3.xls'),
  sheet = 1,
  col_names = TRUE
) %>%
  filter(`dN/dS` > 1) %>%
  select(`Associated Gene Name`) %>%
  pull()

# Peng et al. 2020
peng_471 <- read_csv(
  file = here(
    'data', 'literature-marine-genes', 
    'peng-2020', 'tabula-rough-export-from-pdf.csv'
  ), 
  col_names = c('custom_id', 'symbol', 'description'),
  col_types = cols(), 
  na = ""
) %>%
  select(symbol) %>%
  filter(
    !is.na(symbol),
    nchar(symbol) > 2,
    ! symbol == '0278.1|'
  ) %>%
  mutate(
    symbol = str_replace(
      string = symbol, 
      pattern = 'gi\\|602626649\\|ref\\|XP_00742', 
      replacement = 'gi\\|602626649\\|ref\\|XP_007420278.1'
    ),
    symbol = str_replace(
      string = symbol,
      pattern = 'ETE71823\\.',
      replacement = 'ETE71823.1'
    ),
    symbol = sub('.*\\|', '', symbol),
    symbol = sub('_.*', '', symbol)
  ) %>% 
  pull(symbol)

# Sun et al. 2012 - Dolphin Ensembl IDs do not match Ensembl V66 - can't
# actually convert to a useful identifier...
# sun_693 <- readxl::read_xls(
#   path = here(
#     'data', 'literature-marine-genes','sun-2012',
#     'evs123_Supplementary_Data',
#     'Supplementary_Table_2_PSGs_and_Functinal_clusters.xls'
#   ),
#   sheet = 1,
#   col_names = TRUE
# ) %>%
#   filter(FDR <= 0.05) %>%
#   pull(`Ensembl GeneID`)
#
# dolphin_gtf <- rtracklayer::import(
#   '~/Desktop/Tursiops_truncatus.turTru1.66.gtf.gz', 
# ) %>%
#   as.data.frame() %>%
#   as_tibble() %>%
#   select(gene_id, gene_name)

# Store each dataset as list
lit_genes <- list(
  'chikina' = chikina_690,
  'foote' = foote_191,
  'gayk' = Gayk_152,
  'mcgowan' = mcgowan_228,
  'peng' = peng_471
)

# write_rds(
#   x = lit_genes, 
#   file = 'data/literature-marine-genes/list-literature-genes.rds',
#   compress = 'gz'
# )

# ---------------------------------------------------------------------------- #
# Get gene synonyms from biomaRt
gene_synonyms <- imap_dfr(lit_genes, ~{
  df <- biomaRt::getBM(
    attributes = attributes,
    filters = "hgnc_symbol",
    values = .x,
    mart = ensembl
  ) %>%
    as_tibble()
  
  # Append original query symbols to synonym column - this column is used for searching
  df2 <- 
    tibble(
      hgnc_symbol = .x,
      external_synonym = .x
    )
  
  # Join query tibble with original symbol tibble - filter out junk
  df_bind <- bind_rows(
    df, df2
  ) %>%
    filter(external_synonym != '') %>%
    select(
      symbol = hgnc_symbol,
      synonym = external_synonym
    ) %>%
    arrange(symbol) %>%
    mutate(
      dataset = .y
    )
}) %>%
  distinct() %>%
  mutate(
    symbol = toupper(symbol),
    synonym = toupper(synonym),
    dataset = case_when(
      dataset == 'chikina' ~ 'Chikina et al. 2016',
      dataset == 'foote' ~ 'Foote et al. 2015',
      dataset == 'gayk' ~ 'Gayk et al. 2018',
      dataset == 'mcgowan' ~ 'McGowan et al. 2012',
      dataset == 'peng' ~ 'Peng et al. 2020'
    )
  )

# write_csv(
#   x = gene_synonyms,
#   file = here('data', 'literature-marine-genes', 'literature-marine-genes.csv')
# )

# ---------------------------------------------------------------------------- #
# Join tables to find overlap
fs::dir_create(
  path = here('figures', 'gene-overlap')
)

overlap_379 <- gene_synonyms[gene_synonyms$synonym %in% marine_379_synonyms$synonym, ] %>%
  distinct(symbol, .keep_all = TRUE) %>%
  arrange(symbol)

write_csv(
  x = overlap_379,
  file = here(
    'data', 'selection-results',
    'summary-tables', 'Four-PSGs-overlap-literature.csv'
  ), 
  col_names = TRUE
)

ragg::agg_png(
  filename = here('figures', 'gene-overlap', 'laevis-379-overlap.png'),
  width = 1500,
  height = 1500,
  units = 'px',
  res = 144
)
overlap_379 %>%
  group_by(dataset) %>%
  summarise(count = n()) %>%
  mutate(
    n_genes = case_when(
      dataset == 'Chikina et al. 2016' ~ '690',
      dataset == 'Foote et al. 2015' ~ '191',
      dataset == 'Gayk et al. 2018' ~ '152',
      dataset == 'McGowan et al. 2012' ~ '228',
      dataset == 'Peng et al. 2020'~ '471'
    )
  ) %>%
  ggplot(
    aes(
      x = dataset,
      y = count,
      fill = dataset
    )
  ) +
  geom_bar(stat = 'identity', colour = 'black') +
  geom_text(aes(label = n_genes), vjust = -0.5) +
  ggpomological::scale_fill_pomological() +
  scale_y_continuous(breaks = seq(1,14, 1), limits = c(0, 14)) +
  theme_bw() +
  theme(
    # Axis labels
    axis.title = element_blank(),
    # axis.title.x = element_text(size = 16, face = 'bold'),

    # Axis text
    axis.text = element_text(size = 14),

    # Legend
    # legend.title = element_blank(),
    legend.position = 'none',
    # legend.text = element_text(size = 14)
  )
invisible(dev.off())

# No overlap
gene_synonyms[gene_synonyms$synonym %in% marine_6_synonyms$synonym, ] %>%
  distinct(symbol, .keep_all = TRUE) %>%
  arrange(symbol)
