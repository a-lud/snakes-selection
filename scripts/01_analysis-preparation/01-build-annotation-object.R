# Accessory scripts
source('scripts_paper/util/readFunannotateAnnotations.R')
source('scripts_paper/util/annotateAllOGID.R')
source('scripts_paper/util/annotateGO.R')

# ---------------------------------------------------------------------------- #
# Annotation object for 9 snake dataset
out <- annotateAllOGID(
  path_funannotate = '../../analyses/annotation-station/funannotate', 
  path_omaheaders = 'data/oma-16-samples/Output/OrthologousGroups.csv',
  species_funannotate = c('aipysurusLaevis','najaNaja'), 
  species_ncbi = c('notechisScutatus', 'pseudonajaTextilis'), 
  path_gff3 = '../../analyses/genome-data/annotations',
  path_outfile = 'data/utility-data/master-annotation.rds'
)

# Master Gene Ontology annotation file - GO ids and descriptions relating to orthologue genes
annotateGO(
  annotation = 'data/utility-data/master-annotation.rds', 
  path_outfile = 'data/utility-data/master-GO-annotation.csv'
)


