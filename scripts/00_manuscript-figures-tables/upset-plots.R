# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
  library(ComplexHeatmap)
})

fs::dir_create(
  path = here('manuscript', 'figures')
)

# ---------------------------------------------------------------------------- #
# Import HyPhy data
hyphy <- read_rds(
  file = here('data', 'hyphy-dataframes/hyphy-four-sample.rds')
)

# ---------------------------------------------------------------------------- #
# Significant OMA identifiers belonging to each snake
s1 <- hyphy$absrel$`Branch Attributes` %>%
  filter(Species == 'aipysurusLaevis',
         `Corrected P-value` <= 0.001) %>%
  mutate(ID = sub('.*-', '', ID)) %>%
  pull('ID')

s2 <- hyphy$absrel$`Branch Attributes` %>%
  filter(Species == 'notechisScutatus',
         `Corrected P-value` <= 0.001) %>%
  mutate(ID = sub('.*-', '', ID)) %>%
  extract2('ID')

s3 <- hyphy$absrel$`Branch Attributes` %>%
  filter(Species == 'pseudonajaTextilis',
         `Corrected P-value` <= 0.001) %>%
  mutate(ID = sub('.*-', '', ID)) %>%
  extract2('ID')

s4 <- hyphy$absrel$`Branch Attributes` %>%
  filter(Species == 'najaNaja',
         `Corrected P-value` <= 0.001) %>%
  mutate(ID = sub('.*-', '', ID)) %>%
  extract2('ID')

l <- list(
  'Aipysurus laevis' = s1,
  'Notechis scutatus' = s2,
  'Pseudonaja textilis' = s3,
  'Naja naja' = s4
)

# ---------------------------------------------------------------------------- #
# Combinations matrix + other objects used by UpSet
mat <- make_comb_mat(l)
setsize = set_size(mat)
combsize = comb_size(mat)
combdeg = comb_degree(mat)

# ---------------------------------------------------------------------------- #
# Create the plots using 'UpSet' plot
upset <- UpSet(
  m = mat,
  
  # Data order
  set_order = names(sort(x = setsize, decreasing = TRUE)), # Row order
  comb_order = c(1,2,5,3,4,11,9,6,7,10,8, 14, 13, 15, 12), # Column order
  
  # Annotation settings
  pt_size = unit(6, "mm"), # Size of points representing combinations
  lwd = 2, # Line width
  row_names_gp = grid::gpar(fontface = 'italic', fontsize = 16), # Control row names
  # comb_col = c("#c03728", "#919c4c", "#f5c04a", "#828585")[combdeg], # Colour of lines
  
  # Bar annotations at the top
  top_annotation = HeatmapAnnotation(
    "Intersecting\nGenes\n" = anno_barplot(
      x = combsize,
      ylim = c(0, 380),
      border = FALSE,
      
      # Bar-plot axis parameters
      axis_param = list(
        'at' = seq(0,400,100),
        'gp' = gpar(fontsize = 14)
      ),
      
      # Numbers above columns
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 14),
      
      # Column colours
      # gp = gpar(
      #   fill = c("#c03728", "#919c4c", "#f5c04a", "#828585")[combdeg],
      #   col = c("#c03728", "#919c4c", "#f5c04a", "#828585")[combdeg]
      # ),
      height = unit(4, 'cm')
    ),
    
    # Controlling the bar-plot axis label
    annotation_name_side = 'left',
    annotation_name_rot = 90,
    annotation_name_gp = gpar(
      fontsize = 16,
      fontface = 'bold'
    )
  ),
  
  # Bar plot annotation for set information (no. PSGs per sample)
  right_annotation = rowAnnotation(
    annotation_name_gp = gpar(
      fontsize = 16,
      fontface = 'bold'
    ),
    "\nPSGs Per\nSample" = anno_barplot(
      x = setsize, 
      baseline = 0, 
      which = 'row', 
      border = FALSE,
      # Numbers above columns
      add_numbers = TRUE,
      numbers_gp = gpar(
        fontsize = 14
      ),
      gp = gpar(
        fill = c('#118ab2', '#2a9d8f', '#ffb703', '#e76f51')
        # col = c('#118ab2', '#2a9d8f', '#ffb703', '#e76f51')
      ),
      ylim = c(0, 550), 
      axis_param = list(
        'at' = seq(0,600,100),
        'gp' = gpar(fontsize = 14)
      ),
      width = unit(4, 'cm')
    )
  )
)

# ---------------------------------------------------------------------------- #
# Save plot
ragg::agg_png(
  filename = here('manuscript', 'figures', 'supplementary-upset.png'),
  width = 1000,
  height = 1000,
  units = 'px',
  res = 120
)
print(upset)
invisible(dev.off())
