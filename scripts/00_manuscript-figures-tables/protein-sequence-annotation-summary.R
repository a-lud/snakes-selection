# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
})

fs::dir_create(
  path = here('manuscript', 'tables')
)

# ---------------------------------------------------------------------------- #
# Parse Funannotate annotation tables
funannotate <- fs::dir_ls(
  path = here('..', '..', 'analyses', 'annotation-station', 'funannotate'), 
  recurse = TRUE,
  glob = '*.annotations.txt'
) %>%
  extract(c(1,4)) %>% 
  as.character() %>%
  set_names(sub('.annotations.txt', '', basename(.))) %>%
  map(vroom::vroom, col_types = cols(), col_names = TRUE)

# ---------------------------------------------------------------------------- #
# Import protein sequence summary statistics
sequences <- fs::dir_ls(
  path = here('data', 'utility-data'), 
  glob = '*.stats'
) %>%
  as.character() %>%
  read_tsv(id = 'dataset', show_col_types = FALSE) %>%
  mutate(
    dataset = sub('.stats', '', basename(dataset)),
    dataset = case_when(
      dataset == 'complete-orf' ~ 'Complete Proteins',
      dataset == 'raw-cds' ~ 'Total Proteins'
    ),
    file = sub('\\..*', '', file)
  ) %>%
  select(dataset, file, num_seqs) %>%
  mutate(
    file = case_when(
      file == 'Atenuis_CAGATC_TTAGGC_GCCAAT' ~ 'Aipysurus tenuis',
      file == 'Nscutatus_L1_eye' ~ 'Notechis scutatus (eye)',
      file == 'Ptextilis_L1_eye' ~ 'Pseudonaja textilis (Eye)',
      file == 'aipysurusLaevis' ~ 'Aipysurus laevis',
      file == 'aipysurusLaevis_ACAGTG_ACTTGA_CGATGT' ~ 'Aipysurus laevis (Tail/Body)',
      file == 'aipysurusLaevis_KLS0468' ~ 'Aipysurus laevis (Vomeronasal)',
      file == 'aipysurusLaevis_L1_eye' ~ 'Aipysurus laevis (Eye)',
      file == 'aipysurusMosaicus' ~ 'Aipysurus mosaicus',
      file == 'hydrelapsDarwiniensis' ~ 'Hydrelaps darwiniensis',
      file == 'hydrophisKingii' ~ 'Hydrophis kingii',
      file == 'hydrophisMajor' ~ 'Hydrophis major (Eye)',
      file == 'hydrophisMajor_KLS0460' ~ 'Hydrophis major (Liver/Heart/Testis)',
      file == 'hydrophisPeronii' ~ 'Hydrophis peronii',
      file == 'najaNaja' ~ 'Naja naja',
      file == 'notechisScutatus' ~ 'Notechis scutatus',
      file == 'pseudonajaTextilis' ~ 'Pseudonaja textilis',
      file == 'Alaevis_ACAGTG_ACTTGA_CGATGT' ~ 'Aipysurus laevis (Tail/Body)',
      file == 'Alaevis_KLS0468' ~ 'Aipysurus laevis (Vomeronasal)',
      file == 'Alaevis_L1_eye' ~ 'Aipysurus laevis (Eye)',
      file == 'AMO_L1' ~ 'Aipysurus mosaicus',
      file == 'APE_L1' ~ 'Hydrophis peronii',
      file == 'DKI_L2' ~ 'Hydrophis kingii',
      file == 'HDA_L1' ~ 'Hydrelaps darwiniensis',
      file == 'Hmajor_KLS0460' ~ 'Hydrophis major (Liver/Heart/Testis)',
      file == 'HMA_L2' ~ 'Hydrophis major (Eye)'
    )
  ) %>%
  pivot_wider(names_from = dataset, values_from = num_seqs) %>%
  mutate(
    `Data type` = ifelse(
      file %in% c('Aipysurus laevis', 'Notechis scutatus', 'Pseudonaja textilis', 'Naja naja'),
      'Genome',
      'Transcriptome'
      ),
    `Total Proteins` = ifelse(`Data type` == 'Transcriptome', NA_integer_, `Total Proteins`)
  ) %>%
  arrange(`Data type`, file)

# ---------------------------------------------------------------------------- #
# Number of genes with any functional annotation at all
funannotate.filtered <- imap_dfr(funannotate, ~{
  filtered <- .x %>%
    select(GeneID, Name, PFAM, InterPro, `GO Terms`, Secreted, Membrane, Protease, CAZyme) %>%
    filter_at(
      .vars = vars(2:9),
      .vars_predicate = any_vars(!is.na(.))
    )
  
  total_gene_models <- .x %>% pull(GeneID) %>% unique() %>% length()
  with_gene_symbols <- filtered %>% filter(!is.na(Name)) %>% pull(GeneID) %>% unique() %>% length()
  with_func_annotation <- filtered %>% pull(GeneID) %>% unique() %>% length()
  
  # Summary table
  tribble(
    ~Sample, ~`Predicted Gene Models`, ~`Gene Models with Symbol`, ~`Gene Models With Any Annotation`,
    .y, total_gene_models, with_gene_symbols, with_func_annotation
  )
}) %>%
  mutate(
    Sample = sub('_', ' ', Sample)
  )

# ---------------------------------------------------------------------------- #
# Write gene model and protein summary tables
funannotate.filtered %>%
  write_csv(
    file = here('manuscript', 'tables', 'supplementary-table-funannotate.csv'), 
    col_names = TRUE
  )

sequences %>%
  write_csv(
    file = here('manuscript', 'tables', 'supplementary-table-protein-sequence.csv'),
    col_names = TRUE
  )
