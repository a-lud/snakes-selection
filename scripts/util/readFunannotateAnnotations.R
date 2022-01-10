# -------------------------------------------------------- #
# Function to import '.annotation.txt' files
# 
# Simple function to import Funannotate annotation
# files for a specific subset of annotation columns
#
# -------------------------------------------------------- #
readFunannotateAnnotations <- function(path, species=NULL) {
  
  # Annotation files + subset
  f <- list.files(path = path, 
                  pattern = '.annotations.txt', 
                  full.names = TRUE, 
                  recursive = TRUE)
  
  # Set species names to follow format camel case
  # NOTE: Funannotate annotation file will always follow format of:
  #       species_name.annotation.txt - underscore always between species' name components
  nm <- stringr::str_split(string = sub('.annotations.txt', '', basename(f)), pattern = '_')
  nm <- purrr::map(nm, function(sp){
    sp[1] <- tolower(x = sp[1])
    sp[2] <- stringr::str_to_title(string = sp[2])
    sp <- paste(sp, collapse = '')
  })
  names(f) <- nm
  
  # Subset for species of interest
  if(!is.null(species)) {
    f <- f[names(f) %in% species] # subset for 'species'
  }

  # Columns of interest from Funannotate file
  col_subset <- c('GeneID','Name','PFAM',
                  'InterPro','EggNog','GO Terms')

  # Import data using vroom - fast
  df <- purrr::map(f, ~{
    vroom::vroom(file = .x,
                 delim = '\t',
                 col_names = TRUE,
                 col_types = vroom::cols(),
                 col_select = dplyr::all_of(col_subset))
  })

  # Return list of tibbles - named by species
  return(df)
}
