# Libraries ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
})

# Data ---
an <- read_rds(file = 'data/utility-data/master-annotation.rds')
og <- fs::dir_ls(
  path = 'data/oma-16-samples/complete-oma-groups/four-sample/nucleotide',
  all = TRUE
) %>%
  sub(
    '.*/(.*).fa', '\\1', .
  )

# Only the 3,594 genes that have orthologues (D4)
sub <- an[an$Group %in% og, ]

# Genes of interest - 379 PSGs from A. laevis - removed from candidature
d4 <- read_lines(file = 'data/selection-results/four/absrel-pval-0.001-marine.txt')
sub_no_marine <- sub[! sub$Group %in% d4, ] %>%
  ungroup() %>%
  select(Group, Gene) %>%
  distinct()

# Marine specific genes from the literature
lit_genes <- read_csv(
  file = here('data', 'literature-marine-genes', 'literature-marine-genes.csv'),
  col_names = TRUE,
  col_types = cols()
)

# HyPhy data
hyphy <- read_rds(
  file = here('data', 'hyphy-dataframes', 'hyphy-four-sample.rds')
)

# Get gene symbol synonyms ---
ensembl <- biomaRt::useEnsembl(
  biomart = 'genes',
  dataset = 'hsapiens_gene_ensembl'
)

# Attributes to return
attributes <- c('external_synonym', "hgnc_symbol")

# Query to BiomaRt
synonyms <- biomaRt::getBM(
  attributes = attributes,
  filters = "hgnc_symbol",
  values = sub_no_marine$Gene,
  mart = ensembl
) %>%
  as_tibble()

synonyms <- bind_rows(
  tibble(
    hgnc_symbol = sub_no_marine$Gene,
    external_synonym = sub_no_marine$Gene
  ),
  synonyms
) %>%
  select(
    Gene = hgnc_symbol,
    Synonym = external_synonym
  ) %>%
  left_join(sub_no_marine) %>%
  arrange(Gene) %>%
  filter(
    !str_detect(string = Synonym, pattern = 'C.*orf.*'),
    !str_detect(Synonym, 'LOC.*')
  )

# Intersect with literature genes ---
overlap <- inner_join(lit_genes, synonyms, by = c('synonym' = 'Synonym')) %>% # , 'symbol' = 'Gene'
  distinct(symbol, dataset, .keep_all = TRUE) %>%
  arrange(symbol)

overlap %>%
  ggplot(
    aes(
      x = dataset,
      fill = dataset
    )
  ) +
  geom_histogram(
    stat = 'count',
    colour = 'black'
  ) +
  ggpomological::scale_fill_pomological() +
  labs(
    x = 'Study'
  ) +
  # scale_y_continuous(breaks = seq(1,14, 1), limits = c(0, 14)) +
  theme_bw() +
  theme(
    # Axis labels
    axis.title = element_blank(),
    # axis.title.x = element_text(size = 16, face = 'bold'),
    
    # Axis text
    axis.text = element_text(size = 14),
    
    # Legend
    legend.title = element_blank(),
    legend.position = 'bottom',
    legend.text = element_text(size = 14)
  )

# Get aBSREL and RELAX information for these genes
absrel <- hyphy$absrel$`Branch Attributes` %>%
  separate(col = ID, into = c('model', 'Group'), sep = '-') %>%
  inner_join(overlap) %>%
  # filter(Species == 'aipysurusLaevis') %>%
  select(
    Group, Species, `aBSREL Corrected P-value` = `Corrected P-value`
  ) %>%
  filter(!is.na(`aBSREL Corrected P-value`))

relax <- hyphy$relax$`Test Results` %>%
  separate(col = ID, into = c('model', 'Group'), sep = '-') %>%
  inner_join(overlap) %>%
  filter(
    ! Group %in% d4
  ) %>%
  select(
    Group,
    `RELAX P-value` = `p-value`,
    `Selection Parameter (k)` = `relaxation or intensification parameter`,
    symbol, Gene, dataset
  )

# Join aBSREL + RELAX data
summary_tbl <- full_join(absrel, relax) %>%
  arrange(Group, `aBSREL Corrected P-value`) %>%
  filter(`aBSREL Corrected P-value` <= 0.05) %>%
  mutate(
    `RELAX P-value` = ifelse(
      Species == 'aipysurusLaevis',
      `RELAX P-value`, 
      0
    ),
    `Selection Parameter (k)` = ifelse(
      Species == 'aipysurusLaevis',
      `Selection Parameter (k)`,
      0
    ),
    Species = case_when(
      Species == 'aipysurusLaevis' ~ 'Aipysurus laevis',
      Species == 'najaNaja' ~ 'Naja naja',
      Species == 'pseudonajaTextilis' ~ 'Pseudonaja textilis',
      Species == 'notechisScutatus' ~ 'Notechis scutatus'
    )
  )

summary_tbl %>%
  write_csv(
    file = here('data', 'selection-results', 'alternate-PSGs-absrel-relax.csv'),
    col_names = TRUE
  )

# Plot significant aBSREL results ---
ragg::agg_png(
  filename = here('figures', 'gene-overlap', 'non-PSG-absrel-relax-pval.png'),
  width = 2000,
  height = 2000,
  units = 'px',
  res = 144
)
summary_tbl %>%
  select(Group, Species, `aBSREL Corrected P-value`, symbol) %>%
  distinct() %>%
  ggplot(
    aes(
      x = symbol,
      y = `aBSREL Corrected P-value`,
      colour = Species
    )
  ) + 
  geom_point(size = 3) +
  scale_y_reverse() +
  ggpomological::scale_colour_pomological() +
  labs(
    x = 'Gene'
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size=3)
    )
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14, angle = 90),
    axis.text.y = element_text(size = 14),
    
    # Legend
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = 'italic'),
    legend.position = 'bottom'
  )
invisible(dev.off())
