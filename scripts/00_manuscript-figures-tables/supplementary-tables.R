# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
})

# ---------------------------------------------------------------------------- #
# Output directory
fs::dir_create(
  path = here('manuscript', 'tables'), 
  recurse = TRUE
)

# ---------------------------------------------------------------------------- #
# Import HyPhy data
hyphy <- read_rds(
  file = here('data', 'hyphy-dataframes', 'hyphy-four-sample.rds')
)

# ---------------------------------------------------------------------------- #
# PSG lists
marine <- read_lines(here('data', 'selection-results', 'four', 'absrel-pval-0.001-marine.txt'))
terrestrial <- read_lines(here('data', 'selection-results', 'four', 'absrel-pval-0.001-terrestrial.txt'))

# ---------------------------------------------------------------------------- #
# Annotation
annotation <- read_csv(
  file = here('data', 'utility-data', 'master-GO-annotation.csv.gz'),
  col_names = TRUE, 
  col_types = cols()
) %>%
  select(
    `OMA ID` = Group,
    Gene,
    GO,
    TERM,
    DEFINITION,
    ONTOLOGY
  ) %>%
  distinct()
  # filter(
  #   if_all(
  #     .cols = 2:6,
  #     .fns = ~!is.na(.x)
  #   )
  # )

gene_names <- annotation %>% select(`OMA ID`, Gene) %>% distinct()

# ---------------------------------------------------------------------------- #
# Discussion genes data
genes.list <- list(
  "Respiratory" = c("HIKESHI", "PTGES3", "FOXA1", "GRN", "YWHAZ", "ITGB6", "MECP2"),
  "Circulatory" = c("ARSJ", "MECP2" , "P2RX4", "ADIPOQ", "PRDM1", "C1GALT1", "HAND2", 
                    "YWHAZ", "RANGRF", "BBS7", "FGB", "CCDC137", "LMOD2", "ALAD", 
                    "CYP11B", "LYN", "CIR1", "NCOA4"),
  "Hypoxia Tolerance" = c("ANPEP", "ATG7", "VCPIP1", "ATF4", "GTSE1", "RAD17", "ERCC1", "RPA2", "ASF1A",
                          "RBBP8", "PRDX6", "GPX8", "GSKIP", "EXO1"),
  "Salt and Water Balance" = c("SLC2A1", "SLC16A5", "SLC39A13", "SLC2A3", "TSC1", "C1GALT1", "SIX1A"),
  "Neural and Sensory" = c('CRX', 'PRDM1', 'CNTF', 'CABP4', 'PRDM1', 'ATF4', 'CNTF', 'CABP4', 'PPT1', 'B4GALT2', 'MECP2')
)

# Temporary table: Adaptation, Gene, OMA ID
genes.discussion <-  genes.list %>%
  tibble(Adaptation = names(.), Gene = .) %>%
  unnest(cols = Gene) %>%
  left_join(gene_names %>% filter(`OMA ID` != 'OG48040')) %>%
  distinct(Adaptation, Gene, .keep_all = TRUE)

# ---------------------------------------------------------------------------- #
# GO-DAG summary
gosummaries <- read_rds(
  file = here('data', 'utility-data', 'updatedGOSummary.rds')
) %>%
  filter(
    terminal_node == TRUE,
    shortest_path >= 3,
    longest_path >= 3
  )

# ---------------------------------------------------------------------------- #
# aBSREL data - cleaned
absrel <- hyphy[['absrel']][['Branch Attributes']] %>%
  left_join(hyphy$absrel$Input) %>%
  mutate(ID = sub('.*-', '', ID)) %>%
  unnest(cols = `Rate Distribution`) %>%
  group_by(ID, Species) %>%
  mutate(
    value = paste0(
      c('omega', 'proportion'),
      rep(seq(1:(n()/2)), each = 2)
    )
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = value, values_from = `Rate Distribution`) %>%
  left_join(gene_names, by = c('ID' = 'OMA ID')) %>%
  select(ID, Gene, everything(),
         -c(
           'Partition', 'File Name', 'Number of Sequences',
           'Original Name', 'Partition Count', 'Trees'
         )
  ) %>%
  rename(
    `OMA ID` = ID
  ) %>%
  filter(Species != 'Node2')

# aBSREL subsets: PSG (marine + terrestrial)
absrel.psg.marine <- absrel[absrel$`OMA ID` %in% marine, ] %>%
  arrange(`OMA ID`, `Corrected P-value`)

absrel.psg.terrestrial <- absrel[absrel$`OMA ID` %in% terrestrial, ] %>%
  arrange(`OMA ID`, `Corrected P-value`)

# absrel subsets: Curated (marine only)
absrel.psg.curated.marine <- absrel.psg.marine[
  absrel.psg.marine$`OMA ID`%in% unique(annotation[annotation$GO %in% gosummaries$id, ]$`OMA ID`), 
] %>%
  filter(
    Species == 'aipysurusLaevis'
  ) %>% 
  filter(!str_detect(string = Gene, pattern = '^OG'))

# Write the curated 209 to file (needed for sanity check when manually curating)
write_csv(
  x = absrel.psg.curated.marine, 
  file = here('manuscript', 'tables', 'curated-209-orthologs.csv'), 
  col_names = TRUE
)

# ---------------------------------------------------------------------------- #
# RELAX data - cleaned
relax <- hyphy[['relax']][['Test Results']] %>%
  left_join(hyphy$relax$Input) %>%
  rename(
    K = `relaxation or intensification parameter`,
    `P-value` = `p-value`
  ) %>%
  mutate(
    Condition = ifelse(
      K > 1,
      'Intensification',
      'Relaxation'
    ),
    ID = sub('.*-', '', ID)
  ) %>%
  rename(`OMA ID` = ID) %>%
  left_join(gene_names) %>%
  select(`OMA ID`, Gene, everything(), -c('File Name', 'Number of Sequences', 'Partition Count', 'Trees'))

relax.significant <- relax %>%
  filter(`P-value` <= 0.001)

# ---------------------------------------------------------------------------- #
# Discussion Genes (48) table
discussion.genes.48 <- genes.discussion %>%
  left_join(absrel %>% filter(Species == 'aipysurusLaevis')) %>%
  left_join(relax, by = c('OMA ID', 'Gene')) %>%
  select(
    Adaptation, Gene, 
    `Corrected P-value`, 
    `Full adaptive model (non-synonymous subs/site)`,
    `Full adaptive model (synonymous subs/site)`,
    LRT.x,
    omega1, proportion1, omega2, proportion2,
    K, `P-value`, LRT.y, Condition
  ) %>%
  mutate(
    across(
      .cols = c(4:11, 13),
      .fns = ~round(.x, 3)
    )
  ) %>%
  arrange(Adaptation, `Corrected P-value`)

# ---------------------------------------------------------------------------- #
# Objects to write to file
# - aBSREL - all genes
# - RELAX - all genes
# - Subsets:
#   - absrel PSG
#     - marine [x]
#     - terrestrial [x]
#   - absrel PSG + curated marine [x]
#   - relax significant [x]
# - 48 discussion genes
# - Annotation - all genes

sheets <- list(
  'aBSREL (3,594)' = absrel,
  'aBSREL PSG A. laevis (379)' = absrel.psg.marine,
  'aBSREL PSG Terrestrial (514)' = absrel.psg.terrestrial,
  'aBSREL Curated (209)' = absrel.psg.curated.marine,
  'RELAX (3,543)' = relax,
  'RELAX Significant (418)' = relax.significant,
  'Genes in discussion (48)' = discussion.genes.48,
  'Annotation table (3,594)' = annotation
)

# ---------------------------------------------------------------------------- #
# Write XLSX file with each list element as a sheet
openxlsx::write.xlsx(
  x = sheets, 
  file = here(
    'manuscript', 'tables',
    'supplementary-file-tables.xlsx'
  ), 
  overwrite = TRUE
)


