# -------------------------------------------------------- #
# Create a master annotation file mapping OMA ids to gene
#
# Map annotations to all orthologues. Not filtering on
# Specific samples.
#
# -------------------------------------------------------- #

source('scripts_paper/util/readFunannotateAnnotations.R')

annotateAllOGID <- function(path_funannotate,
                            path_omaheaders,
                            path_gff3 = NULL,
                            species_funannotate=NULL,
                            species_ncbi=NULL,
                            path_outfile = NULL) {
  
  # Import annotation + header files
  lst_annotation <- readFunannotateAnnotations(path = path_funannotate, 
                                               species = species_funannotate)
  df_header <- vroom::vroom(file = path_omaheaders,
                            delim = ',',
                            col_names = TRUE,
                            col_types = vroom::cols())
  
  # Iterate over each species, building annotation table
  df_spec_anno <- purrr::map(colnames(df_header)[-1], function(species){
    
    # Build annotation table based on annotation method
    if (species %in% species_funannotate) {
      # Get species annotation table
      df_anno <- lst_annotation[[species]]
      df_anno <- dplyr::select(df_anno, -EggNog, -PFAM, -InterPro)
      df_anno <- dplyr::rename(df_anno, Gene = Name, GO = `GO Terms`)
      df_anno <- dplyr::mutate(.data = df_anno,
                               GO = stringr::str_split(string = GO, pattern = ';'))
      
      # Clean up species dataframe
      df_spec <- dplyr::select(.data = df_header, Group, dplyr::all_of(species))
      df_spec <- dplyr::rename(.data = df_spec,
                               GeneID = species)
      df_spec <- dplyr::mutate(.data = df_spec, GeneID = sub('-T.*', '', GeneID))
      df_spec <- tidyr::drop_na(data = df_spec)
      
      # # Join Annotation with OGID
      df <- dplyr::left_join(x = df_spec, y = df_anno, by = 'GeneID')
      df <- dplyr::distinct(.data = df, Group, .keep_all = TRUE) 
      df <- dplyr::mutate(.data = df,
                          Gene = ifelse(is.na(Gene), Group, Gene))
      df <- tidyr::unnest(df, cols = 'GO')
      
    } else if (species %in% species_ncbi) {
      
      # Build biomart dataset name from species - ASSUMED <species>_gene_ensembl
      l1 <- substr(species, 1, 1)
      bm_dataset <- paste0(l1, tolower(x = sub("^.*?([A-Z])", "\\1", species)), '_gene_ensembl')
      
      # Clean up species dataframe
      df_spec <- dplyr::select(.data = df_header, Group, dplyr::all_of(species))
      df_spec <- dplyr::rename(.data = df_spec, Name = species)
      df_spec <- dplyr::mutate(.data = df_spec, Name = sub('^rna-', '', Name))
      
      # Import GFF3 files for NCBI species - get gene list
      gff3 <- list.files(path = path_gff3,
                         pattern = '.gff3',
                         full.names = TRUE,
                         recursive = TRUE)
      names(gff3) <- sub('.gff3', '', basename(gff3))
      
      # Filter for species of interest (basename of gff3)
      gff3 <- gff3[ names(gff3) %in% species ]
      
      # Read GFF3 data
      print('Reading GFF3')
      gff3 <- rtracklayer::readGFF(filepath = gff3,
                                   version = 3,
                                   columns = c('seqid', 'type', 'attributes'))
      
      # Get gene + mrna ids
      cols_gff <- c('type', 'Name', 'gene')
      
      df_gff <- tibble::as_tibble(x = gff3)
      df_gff <- dplyr::select(.data = df_gff, dplyr::all_of(cols_gff))
      df_gff <- dplyr::filter(.data = df_gff, type == 'mRNA')
      df_gff <- dplyr::left_join(x = df_spec, y = df_gff, by = 'Name')
      
      genes <- unique(df_gff[['gene']])
      
      mart <- biomaRt::useMart(biomart = 'ensembl', dataset = bm_dataset)
      print('Mart got')
      
      res <- tibble::as_tibble(biomaRt::getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id', 
                                                             'entrezgene_id', 'go_id', 'go_linkage_type'),
                                              filters = c('hgnc_symbol'),
                                              values = genes,
                                              mart = mart))
      print('BM got')
      
      res <- dplyr::full_join(x = df_gff,
                              y = res,
                              by = c('gene' = 'hgnc_symbol'))
      
      res <- dplyr::select(.data = res, -type, Name)
      res <- dplyr::rename(.data = res,
                           GeneID = Name,
                           Gene = gene,
                           Ensembl = ensembl_gene_id,
                           Entrez = entrezgene_id,
                           GO = go_id,
                           Evidence = go_linkage_type)
    }
  })
  
  # Get unique list entries for PFAM, InterPro, GO
  names(df_spec_anno) <- colnames(df_header)[-1]
  df_spec_anno <- dplyr::bind_rows(df_spec_anno)
  df_spec_anno <- dplyr::group_by(.data = df_spec_anno, Group)
  df_spec_anno <- dplyr::distinct(df_spec_anno, Group, Gene, GO, Evidence, .keep_all = TRUE)
  df_spec_anno <- dplyr::ungroup(x = df_spec_anno)

  # Standardise gene-symbols
  df_anno <- dplyr::mutate(.data = df_spec_anno,
                           Gene = sub('_\\d', '', Gene))
  df_anno <- dplyr::group_by(.data = df_anno, Group)
  df_anno <- dplyr::mutate(.data = df_anno,
                           many_symbol = list(unique(as.character(na.omit(Gene)))),
                           many_symbol = ifelse(length(many_symbol) > 1,
                                                list(unique(unlist(many_symbol)[unlist(many_symbol) != unique(Group)])),
                                                many_symbol),
                           many_symbol = ifelse(length(many_symbol) > 1,
                                                list(unique(unlist(many_symbol)[!stringr::str_detect(unlist(many_symbol), 'LOC\\d{9}')])),
                                                many_symbol),
                           many_symbol = ifelse(rlang::is_empty(unlist(many_symbol)),
                                                list(unique(unlist(Group))),
                                                many_symbol),
                           many_symbol = ifelse(length(many_symbol) > 1,
                                                list(sort(unlist(many_symbol))[1]),
                                                many_symbol),
                           many_symbol = unlist(many_symbol))
  df_anno <- dplyr::select(.data = df_anno,
                           Group, -Gene, many_symbol, Ensembl, Entrez, GO, Evidence)
  df_anno <- dplyr::rename(.data = df_anno, Gene = many_symbol)
  
  # Write to file
  if(!is.null(path_outfile)) {
    readr::write_rds(x = df_anno,
                     file = path_outfile,
                     compress = 'gz')
  }
  
  return(df_anno)
}
