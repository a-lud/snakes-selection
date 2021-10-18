# -------------------------------------------------------- #
# Annotate GO Terms with simple names and descriptions
#
# Using the `GO.db` database, get annotations for GO terms
#
# -------------------------------------------------------- #

annotateGO <- function(annotation, path_outfile = NULL) {
  if(is.character(annotation)) {
    df_anno <- readr::read_rds(file = annotation)
  } else if(is(annotation, 'tbl_df')) {
    df_anno <- annotation
  }
  
  df_go <- AnnotationDbi::select(GO.db::GO.db,
                               columns = c("GOID",
                                           "TERM",
                                           "DEFINITION",
                                           "ONTOLOGY"),
                               keys = unique(as.character(na.omit(df_anno[['GO']]))))
  df_go <- tibble::as_tibble(x = df_go)

  df <- dplyr::left_join(df_anno, df_go,
                         by = c('GO' = 'GOID'))

  readr::write_csv(x = df,
                   file = path_outfile, na = '',
                   col_names = TRUE)
}
