# ---------------------------------------------------------------------------- #
# Over-representation Analysis: topGO
#
# This is a script to test for over-representation of GO terms within the D4
# dataset of genes with a signature of selection. The software topGO is used
# as it implements a elmintation/weighted algorithm for dealing with the
# relatedness of GO terms in the GO-DAG. Selection tests are only applied
# to the A. laevis genes as there are enough under selection to actually
# test for enrichment.
#
#                                 NOT USED ---
# The Terrestrial condition has many genes shared over all three snakes, however
# very few genes are actually shared by all three snakes (the bulk of the genes
# under selection belong to N. naja)
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
  library(topGO)
})

# ---------------------------------------------------------------------------- #
# Data import
go <- read_csv(
  file = here('data', 'utility-data', 'master-GO-annotation.csv.gz'),
  col_names = TRUE,
  col_types = cols()
)

# All orthologues in the four-snake condition
four_sample_genes <- fs::dir_ls(
  path = here(
    'data', 'oma-16-samples',
    'complete-oma-groups', 'four-sample', 'nucleotide'
  )
) %>% 
  basename() %>%
  sub('.fa', '', .)

# Marine specific genes
marine <- read_lines(
  here('data', 'selection-results', 'four', 'absrel-pval-0.001-marine.txt')
)

# ---------------------------------------------------------------------------- #
# Write topGO map file 
# <gene>\t<GOID, GOID, GOID, ..., GOID>\n
go_filtered <- go %>%
  filter(
    !is.na(GO),
    Group %in% four_sample_genes,
  ) %>%
  distinct() %>%
  mutate(ogid_gene = paste(Group, Gene, sep = ':'))

# Separate step as I build different object below
go_filtered %>%
  group_by(ogid_gene) %>%
  summarise(GO = paste(unique(GO), collapse = ', ')) %>%
  distinct(ogid_gene, .keep_all = TRUE) %>%
  mutate(
    GO = sub('^, ', '', GO),
    GO = ifelse(GO == '', NA_character_, GO)
  ) %>%
  filter(!is.na(GO)) %>%
  readr::write_tsv(
    file = here('data', 'utility-data', 'topGo-four-snake-ogidGene-map.tsv'),
    col_names = FALSE
  )

# ---------------------------------------------------------------------------- #
# Form matching vector for marine genes as GO-map
marine <- go_filtered %>%
  filter(
    Group %in% marine
  ) %>%
  pull(ogid_gene) %>%
  unique()

# Don't use either of the below - but can if you want
# gomap_ogid <- go_filtered %>%
#   group_by(Group) %>%
#   summarise(GO = paste(unique(GO), collapse = ', ')) %>%
#   distinct(Group, .keep_all = TRUE) %>%
#   mutate(
#     GO = sub('^, ', '', GO),
#     GO = ifelse(GO == '', NA_character_, GO)
#   ) %>%
#   filter(!is.na(GO)) %>%
#   readr::write_tsv(
#     file = 'data/utility-data/topGo-four-snake-ogid-map.tsv',
#     col_names = FALSE
#   )
#   
# gomap_symbol <- go_filtered %>%
#   group_by(Gene) %>%
#   summarise(GO = paste(unique(GO), collapse = ', ')) %>%
#   distinct(Gene, .keep_all = TRUE) %>%
#   mutate(
#     GO = sub('^, ', '', GO),
#     GO = ifelse(GO == '', NA_character_, GO)
#   ) %>%
#   filter(!is.na(GO)) %>%
#   readr::write_tsv(
#     file = 'data/utility-data/topGo-four-snake-symbol-map.tsv',
#     col_names = FALSE
#   )

# ---------------------------------------------------------------------------- #
# TopGodata Gene list
geneID2GO <- topGO::readMappings(
  here('data', 'utility-data', 'topGo-four-snake-ogidGene-map.tsv')
)
geneNames <- names(geneID2GO)

# ---------------------------------------------------------------------------- #
# A.laevis analysis
fs::dir_create(here('data', 'GO-results'))
path <- here('data', 'GO-results', 'topGO-marine-summary.csv')
if(fs::file_exists(path = path)) {
  fs::file_delete(path = path)
  write(x = c('Ontology,Feasible Genes,Sig. Genes'),  
        file = path, 
  )
} else {
  write(x = c('Ontology,Feasible Genes,Sig. Genes'),  
        file = path, 
  )
}

geneList <- factor(as.integer(geneNames %in% marine))
names(geneList) <- geneNames

# Run a test for each Ontology
topGO_marine <- purrr::map(c('BP', 'MF', 'CC'), ~{
  
  # Create GOData object
  godata <- new(
    "topGOdata", 
    ontology = .x, 
    allGenes = geneList,
    nodeSize = 5,
    annot =  annFUN.gene2GO, 
    gene2GO = geneID2GO
  )
  
  # Write counts to file
  write(
    x = paste(c(.x,
                length(genes(godata)),
                length(sigGenes(godata))),
              collapse = ','),
    file = path,
    append = TRUE
  )
  
  # Weighted algorthim (default - considering hierarchy)
  result_wf <- runTest(
    godata, 
    algorithm = 'weight01', 
    statistic = 'fisher'
  )
  
  # Build summary table
  out <- GenTable(
    godata,
    weightFisher = result_wf,
    orderBy = 'weightFisher',
    topNodes = length(usedGO(godata))
  ) %>%
    tibble::as_tibble() %>%
    mutate(
      adjp = p.adjust(weightFisher, method = 'bonferroni'),
      FDR = p.adjust(weightFisher, method = 'fdr')
    ) %>%
    filter(weightFisher <= 0.05)
  
  # Significant genes (?)
  go_gene_list <- scoresInTerm(
    godata,
    out$GO.ID,
    use.names = TRUE
  )
  
  df_genes <- imap(go_gene_list, ~{
    GO <- .y
    GENES <- names(.x[.x > 1 ]) %>%
      sub('OG.*:', '', .) %>%
      paste(., collapse = ', ') %>%
      unique()
    
    tibble(
      'GO.ID' = GO,
      'Genes' = GENES
    )
  }) %>%
    bind_rows()
  
  out <- left_join(out, df_genes)
  
  return(list(
    'godata' = godata,
    'df' = out
  ))
}) %>%
  set_names(value = c('BP', 'MF', 'CC'))

# ---------------------------------------------------------------------------- #
# Write output
sig_table <- map(topGO_marine, ~{
  .x[['df']]
}) %>%
  bind_rows(.id = 'Ontology') %>%
  left_join(go, by = c('GO.ID' = 'GO')) %>%
  dplyr::select(
    GO = `GO.ID`,
    Ontology,
    Term,
    Annotated, 
    Significant,
    Pvalue = weightFisher,
    `Adj. Pvalue` = adjp,
    FDR,
    Genes,
    Definition = DEFINITION
  ) %>%
  distinct() %>%
  arrange(Ontology, Pvalue)

# Write to file
write_csv(
  x = sig_table,
  file = here('data', 'GO-results', 'topGO-marine.csv'), 
  col_names = TRUE
)

# # ---------------------------------------------------------------------------- #
# # TopGodata prep - Terrestrial
# path <- file.path('data/GO-results', paste0('topGO-terrestrial-topGOdata-summary.csv'))
# if(fs::file_exists(path = path)) {
#   fs::file_delete(path = path)
#   write(x = c('Ontology,Feasible Genes,Sig. Genes'),  
#         file = path, 
#   )
# } else {
#   write(x = c('Ontology,Feasible Genes,Sig. Genes'),  
#         file = path, 
#   )
# }
# geneList <- factor(as.integer(geneNames %in% terrestrial))
# names(geneList) <- geneNames
# 
# # Run a test for each Ontology
# topGO_terrestrial <- purrr::map(.x = c('BP', 'MF', 'CC'), ~{
#   
#   # Create GOData object
#   godata <- new("topGOdata", 
#                 ontology = .x, 
#                 allGenes = geneList,
#                 nodeSize = 5,
#                 annot =  annFUN.gene2GO, 
#                 gene2GO = geneID2GO)
#   
#   # Write counts to file
#   write(
#     x = paste(c(.x, 
#                 length(genes(godata)),
#                 length(sigGenes(godata))), 
#               collapse = ','), 
#     file = path,
#     append = TRUE
#   )
#   
#   # Weighted (considering hierarchy)
#   result_wf <- runTest(godata, algorithm = 'weight01', statistic = 'fisher')
#   
#   # Build summary table
#   out <- GenTable(godata, 
#                   weightFisher = result_wf,
#                   orderBy = 'weightFisher',
#                   topNodes = length(usedGO(godata))) %>%
#     tibble::as_tibble() %>%
#     filter(
#       Annotated >= 5
#     ) %>%
#     mutate(adjp = p.adjust(weightFisher, method = 'bonferroni'),
#            FDR = p.adjust(weightFisher, method = 'fdr')) %>%
#     filter(weightFisher <= 0.05)
#   
#   # Subset annotation object by significant genes used in analysis
#   df_sigGenes_anno <- go %>%
#     filter(Group %in% sigGenes(godata))
#   
#   return(list(
#     'godata' = godata,
#     'df' = out,
#     'genes' = genes(godata),
#     'sigGenes' = sigGenes(godata),
#     'df_sigGenes_anno' = df_sigGenes_anno
#   ))
# }) %>%
#   set_names(value = c('BP', 'MF', 'CC'))
# 
# # ---------------------------------------------------------------------------- #
# # List of topGO results
# topGO_results <- list(
#   'marine' = topGO_marine,
#   'terrestrial' = topGO_terrestrial
# )

# ---------------------------------------------------------------------------- #
# Bind Significance Tables together
# topGO_results_filtered <- imap(topGO_results, ~{
#   map(.x, function(ont) {
#     ont[['df']]
#   }) %>%
#     bind_rows(.id = 'Ontology')
# })
# 
# # Write to file
# iwalk(topGO_results_filtered, ~{
#   fs::dir_create(path = 'data/GO-results')
#   path <- file.path('data/GO-results', paste0('topGO-', .y, '.csv'))
#   write_csv(
#     x = .x,
#     file = path, 
#     col_names = TRUE
#   )
# })
