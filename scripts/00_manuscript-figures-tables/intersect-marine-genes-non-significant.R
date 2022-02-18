# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
})

# ---------------------------------------------------------------------------- #
# Data
an <- read_rds(file = here('data', 'utility-data', 'master-annotation.rds'))
oma.id <- fs::dir_ls(
  path = here('data' , 'oma-16-samples', 'complete-oma-groups', 'four-sample', 'nucleotide'),
  all = TRUE
) %>%
  sub(
    '.*/(.*).fa', '\\1', .
  )

# ---------------------------------------------------------------------------- #
# Only the 3,594 genes that have orthologues
oma.id.subset <- an[an$Group %in% oma.id, ]

# Genes of interest - 379 PSGs from A. laevis - removed from candidature
orthologs <- read_lines(file = here('data', 'selection-results', 'four', 'absrel-pval-0.001-marine.txt'))
sub.no.marine <- oma.id.subset[! oma.id.subset$Group %in% orthologs, ] %>%
  ungroup() %>%
  select(Group, Gene) %>%
  distinct()

# Marine specific genes from the literature
lit.genes <- read_csv(
  file = here('data', 'literature-marine-genes', 'literature-marine-genes.csv'),
  col_names = TRUE,
  col_types = cols()
)

# ---------------------------------------------------------------------------- #
# HyPhy data
hyphy <- read_rds(
  file = here('data', 'hyphy-dataframes', 'hyphy-four-sample.rds')
)

# ---------------------------------------------------------------------------- #
# Get gene symbol synonyms
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
  values = sub.no.marine$Gene,
  mart = ensembl
) %>%
  as_tibble()

# ---------------------------------------------------------------------------- #
# Add in original gene names
synonyms <- bind_rows(
  tibble(
    hgnc_symbol = sub.no.marine$Gene,
    external_synonym = sub.no.marine$Gene
  ),
  synonyms
) %>%
  select(
    Gene = hgnc_symbol,
    Synonym = external_synonym
  ) %>%
  left_join(sub.no.marine) %>%
  arrange(Gene) %>%
  filter(
    !str_detect(string = Synonym, pattern = 'C.*orf.*'),
    !str_detect(Synonym, 'LOC.*')
  )

# ---------------------------------------------------------------------------- #
# Get exact and synonymous matches

# 1. Find exact matches - 141 genes
exact.matches <- synonyms %>%
  select(Gene, Group) %>%
  distinct() %>%
  arrange(Gene) %>%
  inner_join(lit.genes %>% select(symbol, dataset) %>%distinct(), by = c('Gene' = 'symbol'))

# 2. Get matches based on synonym - remove the matches above (as they're accounted for) - and write
#    the file to be checked manually
synonyms.no.exact <- synonyms %>% filter(! Gene %in% exact.matches$Gene)
lit.genes.no.exact <- lit.genes %>% filter( ! lit.genes$symbol %in% exact.matches$Gene )

synonym.matches <- inner_join(
  synonyms.no.exact, lit.genes.no.exact,
  by = c('Synonym' = 'synonym')
) %>% 
  arrange(Gene)

write_csv(
  x = synonym.matches,
  file = here(
    'data', 'literature-marine-genes', 'non-psg-synonym-check.csv'
  ),
  col_names = TRUE
)

# ---------------------------------------------------------------------------- #
# Import genes that matched on synonyms that have been manually checked
synonym.matches.checked <- read_csv(
  file = here('data', 'literature-marine-genes', 'non-psg-synonym-checked.csv'),
  col_names = TRUE,
  col_types = cols()
) %>%
  select(Gene, Group, dataset) %>%
  mutate(dataset = sub('\\\\n', '\n', dataset))

# ---------------------------------------------------------------------------- #
# Bind exact matches with manually curated list
# lit.overlap.non.psg <- bind_rows(exact.matches, synonym.matches.checked)
# lit.overlap.non.psg %>% pull(Gene) %>% unique() %>% length() # 148 Genes
# 
# write_csv(
#   lit.overlap.non.psg,
#   file = here('data', 'literature-marine-genes', 'overlap-non-sig-literature.csv'),
#   col_names = TRUE
# )

lit.overlap.non.psg <- read_csv(here('data', 'literature-marine-genes', 'overlap-non-sig-literature.csv'))

# Match levels to other plot
lit.overlap.non.psg %<>%
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
    dataset = factor(dataset, levels = c(
      'Chikina et al. 2016\n(Marine mammals)',
      'Foote et al. 2015\n(Marine mammals)',
      'Gayk et al. 2018\n(Diving bird)',
      'Peng et al. 2020\n(Sea snake)',
      'McGowan et al. 2012\n(Dolphin)'
    ))
  )

# ---------------------------------------------------------------------------- #
# # Plot
# ragg::agg_png(
#   filename = here('manuscript', 'figures', 'supplementary-literature-overlap-non-psg.png'),
#   width = 1500,
#   height = 800,
#   units = 'px',
#   res = 144
# )
p <- lit.overlap.non.psg %>%
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
  ylim(c(0, 50)) +
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
    
    # Legend
    legend.position = 'none',
  )
write_rds(
  x = p, 
  file = here('data/literature-marine-genes/supplementary-non-psg-lit-overlap.rds')
)
# invisible(dev.off())

# ---------------------------------------------------------------------------- #
# Get aBSREL results for genes
absrel <- hyphy$absrel$`Branch Attributes` %>%
  mutate(
    ID = sub('.*-', '', ID)
  ) %>%
  inner_join(lit.overlap.non.psg, by = c('ID' = 'Group')) %>%
  # filter(Species == 'aipysurusLaevis') %>%
  select(
    ID, Species, `aBSREL LRT` = LRT, `aBSREL Corrected P-value` = `Corrected P-value`, `Full adaptive model`,
    `Full adaptive model (non-synonymous subs/site)`, `Full adaptive model (synonymous subs/site)`,
    `Rate Distribution`, `Rate Classes`, Gene, dataset
  ) %>%
  filter(Species != 'Node2') %>%
  unnest(cols = `Rate Distribution`) %>%
  group_by(ID, Species, Gene, dataset) %>% 
  mutate(
    value = paste0(
      c('omega', 'proportion'),
      rep(seq(1:(n()/2)), each = 2)
    )
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = value, values_from = `Rate Distribution`)

# ---------------------------------------------------------------------------- #
# Get RELAX results for genes
relax <- hyphy$relax$`Test Results` %>%
  mutate(
    ID = sub('.*-', '', ID)
  ) %>%
  inner_join(lit.overlap.non.psg, by = c('ID' = 'Group')) %>%
  select(
    ID, `RELAX LRT` = LRT,
    `RELAX P-value` = `p-value`,
    `Selection Parameter (k)` = `relaxation or intensification parameter`,
    Gene, dataset
  )

# ---------------------------------------------------------------------------- #
# Write xlsx file with results in it
sheets <- list(
  'aBSREL literature overlap' = absrel,
  'RELAX literature overlap' = relax
)

openxlsx::write.xlsx(
  x = sheets, 
  file = here(
    'manuscript', 'tables',
    'supplementary-file-non-psg-literature-overlap.xlsx'
  ), 
  overwrite = TRUE
)
