# ---------------------------------------------------------------------------- #
# RELAX figures + Figure 2

# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
  library(ComplexHeatmap)
  library(patchwork)
})

fs::dir_create(
  path = here('manuscript', 'tables')
)

fs::dir_create(
  path = here('manuscript', 'figures')
)

# ---------------------------------------------------------------------------- #
# Import data
hyphy <- read_rds(
  file = here('data', 'hyphy-dataframes', 'hyphy-four-sample.rds')
)

# Gene lengths
df_length <- hyphy[['relax']][['Input']] %>%
  select(ID, `Number of Sites`)

# PSGs for A. laevis and Terrestrial snakes
marine <- read_lines(
  here('data', 'selection-results', 'four', 'absrel-pval-0.001-marine.txt')
)

terrestrial <- read_lines(
  here('data', 'selection-results', 'four', 'absrel-pval-0.001-terrestrial.txt')
)

# ---------------------------------------------------------------------------- #
# RELAX data

# Format data and filter for significant results (p <= 0.001)
relax_sig <- hyphy[['relax']][['Test Results']] %>%
  left_join(df_length) %>%
  rename(K = `relaxation or intensification parameter`) %>%
  mutate(
    Condition = ifelse(
      K > 1,
      'Intensification',
      'Relaxation'
    ),
    Significance = ifelse(
      `p-value` <= 0.001,
      'Significant',
      'Insignificant'
    ),
    ID = sub('.*-', '', ID),
    # This shows how relax signatures compare to absrel results for both
    # marine and terrestrial gene sets.
    Group = 
      ifelse(
        ID %in% marine,
        'Aipysurus laevis',
        ifelse(
          ID %in% terrestrial,
          'Terrestrial',
          'Remaining'
        )
      )
  ) %>%
  mutate(
    Group = factor(x = Group, levels = c('Aipysurus laevis', 'Terrestrial', 'Remaining')),
    Significance = factor(x = Significance, levels = c('Significant', 'Insignificant'))
  )

# Write data object to file
write_csv(
  x = relax_sig,
  file = here('manuscript', 'tables', 'RELAX-significant-table.csv'), 
  col_names = TRUE
)

# ---------------------------------------------------------------------------- #
# Build the base plot
relax_base <- relax_sig %>%
  ggplot(aes(x = Condition, fill = Condition)) +
  geom_histogram(stat = 'count', colour = 'black', alpha = 0.5) +
  scale_fill_manual(values = c('red', 'black')) +
  labs(
    y = 'Number of Orthologs'
  ) +
  theme_bw() +
  theme(
    # Legend
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = 'none',
    
    # Strip
    strip.text = element_text(size = 16, face = 'bold'),
    
    # Axis title/text
    axis.title = element_text(size = 16, face = 'bold'),
    axis.text = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Significant RELAX results vs insignificant - compare intens/relax (all relative to A. laevis)
ragg::agg_png(
  filename = here('manuscript', 'figures', 'relax-bar-plot-summary-significant-signal.png'),
  width = 1500,
  height = 1000,
  units = 'px',
  res = 144
)
relax_base +
  facet_wrap(
    Significance ~ .,
    scales = 'free_y',
    ncol = 2
  )
invisible(dev.off())

# ---------------------------------------------------------------------------- #
# Facet the significant results by A. laevis vs. Terrestrial
tmp <- relax_sig %>% # Italicise Aipysurus laevis
  filter(Significance == 'Significant')
  
levels(tmp$Group) <- c(
  'Aipysurus laevis' = expression(italic("Aipysurus laevis")),
  'Terrestrial' = 'Terrestrial',
  'Remaining' = 'Remaining'
)

# Faceting controls the overlapping and non-overlapping
p <- tmp %>%
  ggplot(aes(x = Condition, fill = Condition)) +
  geom_histogram(stat = 'count', colour = 'black', alpha = 0.5) +
  scale_fill_manual(values = c('red', 'black')) +
  labs(
    y = 'Number of Orthologs'
  ) +
  theme_bw() +
  theme(
    # Legend
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = 'none',
    
    # Strip
    strip.text = element_text(size = 16, face = 'bold'),
    
    # Axis title/text
    axis.title = element_text(size = 16, face = 'bold'),
    axis.text = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  facet_grid(
    . ~ Group,
    scales = 'free_y', 
    labeller = label_parsed
  )

# ---------------------------------------------------------------------------- #
# UpSet plot of RELAX genes against PSGs belonging to A. laevis and Terrestrial snakes
all_relax <- relax_sig %>% filter(Significance == 'Significant') %>% pull(ID)

l <- list(
  "Significant\nRELAX" = all_relax,
  'Aipysurus laevis\nPSGs' = marine,
  'Terrestrial\nPSGs' = terrestrial
)

mat <- make_comb_mat(l)
setsize <- set_size(mat)
combsize <- comb_size(mat)
combdeg <- comb_degree(mat)

# UpSet can take all the same arguments as complexHeatmap
upset <- UpSet(
  m = mat,
  
  row_labels = gt_render(
    c('Significant<br>RELAX', "*Aipysurus laevis*<br>PSGs", 'Terrestrial<br>PSGs') # Using this to get italics on Aipysurus laevis
  ),
  
  # Data order
  # set_order = names(sort(x = setsize, decreasing = TRUE)), # Row order
  set_order = names(l),
  
  # Annotation settings
  pt_size = unit(6, "mm"), # Size of points representing combinations
  lwd = 2, # Line width
  row_names_gp = grid::gpar(
    # fontface = 'bold',
    fontsize = 16
  ), # Control row names
  
  # Bar annotations at the top
  top_annotation = HeatmapAnnotation(
    "Intersecting\nGenes\n" = anno_barplot(
      x = combsize,
      ylim = c(0, 500),
      border = FALSE,
      
      # Bar-plot axis parameters
      axis_param = list(
        'at' = seq(0,500,100),
        'gp' = gpar(fontsize = 14)
      ),
      
      # Numbers above columns
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 14),
      
      # Column colours
      gp = gpar(
        fill = c('#118ab2', '#919c4c', '#f5c04a', 'grey50', 'grey50')
        # col = c('#118ab2', '#919c4c', '#f5c04a', 'grey50', 'grey50')
      ),
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
  
  # # Bar plot annotation for set information (no. PSGs per sample)
  right_annotation = rowAnnotation(
    annotation_name_gp = gpar(
      fontsize = 16,
      fontface = 'bold'
    ),
    "\nTotal Genes\nPer Set" = anno_barplot(
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
        fill = 'grey50'
      ),
      ylim = c(0, 450),
      axis_param = list(
        'at' = seq(0,500,100),
        'gp' = gpar(fontsize = 14)
      ),
      width = unit(4, 'cm')
    )
  )
)

# Convert to a gg object
ggupset <- ggplotify::as.ggplot(upset)

# ---------------------------------------------------------------------------- #
# Colouring facet labels to match upset
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c('#118ab280', '#919c4c80', '#f5c04a80')
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

# Save plot
ragg::agg_png(
  filename = here('manuscript', 'figures', 'figure-2.png'),
  width = 1500,
  height = 1500,
  units = 'px',
  res = 144
)
ggupset / g + plot_annotation(tag_levels = 'A')
invisible(dev.off())

