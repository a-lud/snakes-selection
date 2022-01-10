# ---------------------------------------------------------------------------- #
# Libraries
options(scipen = 999)
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
  library(patchwork)
  library(ggtree)
})

fs::dir_create(
  path = here('manuscript', 'tables')
)

fs::dir_create(
  path = here('manuscript', 'figures')
)

# ---------------------------------------------------------------------------- #
# Colour palette
colPal <- colorRampPalette(colors = RColorBrewer::brewer.pal(n = 11, name =  'Spectral'))

# ---------------------------------------------------------------------------- #
# Manually curated genes
curated.manual.annotation <- read_csv(
  file = here('manuscript', 'tables', 'curated-209-orthologs-manual-annotations.csv'),
  col_names = TRUE,
  col_types = cols()
) %>%
  select(Gene, System)

max_cols <- max(str_count(string = curated.manual.annotation$System, pattern = '\\/'), na.rm = TRUE) + 1
curated.manual.annotation <- separate(
  data = curated.manual.annotation,
  col = System, sep = '\\/', into = paste0('X', 1:max_cols)
) %>%
  pivot_longer(
    cols = 2:5, names_to = 'blah', values_to = 'systems'
  ) %>%
  filter(
    !is.na(systems)
  ) %>%
  select(-blah) %>%
  mutate(
    broad_grouping = case_when(
      str_detect(systems, 'Cell') ~ 'Cellular functions',
      str_detect(systems, 'Neural|Nervous') ~ 'Neural function/behaviour',
      systems == 'Sensory' ~ 'Sensory perception',
      TRUE ~ as.character(systems)
    )
  ) %>%
  group_by(broad_grouping) %>%
  summarise(
    n_genes = n()
  ) %>%
  arrange(-n_genes) %>%
  slice(1:15) %>%
  mutate(
    col_text = glue::glue("{broad_grouping} ({n_genes})"),
    col_text = factor(x = col_text, levels = col_text),
    broad_grouping = factor(x = broad_grouping, levels = broad_grouping),
    icon = c(
      "<img src='figures/icons/cell.png' width='20' />",
      "<img src='figures/icons/brain.png' width='20' />",
      "<img src='figures/icons/artery.png' width='20' />",
      "<img src='figures/icons/metabolism.png' width='20' />",
      "<img src='figures/icons/immune-system.png' width='20' />",
      "<img src='figures/icons/view.png' width='20' />",
      "<img src='figures/icons/molecules.png' width='20' />",
      "<img src='figures/icons/spine.png' width='20' />",
      "<img src='figures/icons/dna.png' width='20' />",
      "<img src='figures/icons/lungs.png' width='20' />",
      "<img src='figures/icons/salt.png' width='20' />",
      "<img src='figures/icons/muscles.png' width='20' />",
      "<img src='figures/icons/petri-dish.png' width='20' />",
      "<img src='figures/icons/recycle.png' width='20' />",
      "<img src='figures/icons/mitochondria.png'  width='20' />"
    )
  ) %>%
  slice(1:15)

# ---------------------------------------------------------------------------- #
# HyPhy data
hyphy <- read_rds(
  file = here('data', 'hyphy-dataframes', 'hyphy-four-sample.rds')
)

test.branches <- hyphy$absrel$Tested %>% 
  filter(Status == 'test') %>% 
  pull(Node) %>% 
  unique()

# ---------------------------------------------------------------------------- #
# Tree files
tree <- treeio::read.newick(file = here('data', 'trees', 'plotting-d4.nwk')) %>%
  tidytree::as_tibble() %>%
  mutate(
    label = sub('\\{.*', '', label),
    label = case_when(
      label == 'aipysurusLaevis' ~ 'Aipysurus laevis',
      label == 'notechisScutatus' ~ 'Notechis scutatus',
      label == 'pseudonajaTextilis' ~ 'Pseudonaja textilis',
      label == 'najaNaja' ~ 'Naja naja'
    ),
    type = ifelse(is.na(label), NA_character_, 'test'),
    data_type = 'Genome'
  )

# ---------------------------------------------------------------------------- #
# Tree annotation data
comb <- read_csv(
  file = here('manuscript', 'tables', 'combination-matrix.csv'),
  col_names = TRUE,
  col_types = cols()
) %>%
  filter(intersect %in% c(379, 367, 64, 45)) %>%
  pivot_longer(names_to = 'label', values_to = 'bin', -intersect) %>%
  filter(bin != 0) %>%
  select(-bin) %>%
  mutate(
    label = factor(x = label, levels = c(
      'Notechis scutatus',
      'Aipysurus laevis',
      'Pseudonaja textilis',
      'Naja naja'
    )),
    rate = case_when(
      intersect == 367 ~ signif(367/32, digits = 3),
      intersect == 45 ~ signif(45/22, digits = 3),
      intersect == 379 ~ signif(379/18, digits = 3),
      intersect == 64 ~ signif(64/18, digits = 3)
    )
  )

# ---------------------------------------------------------------------------- #
# Join PSG counts to tree data
tree.data <- full_join(tree, comb)

# Convert back to tree object
tree.data <- tidytree::as.treedata(tree = tree.data)

# ---------------------------------------------------------------------------- #
# Horizontal barplot to join with tree
barplot.annotation <- curated.manual.annotation %>% 
  ggplot(
    aes(
      x = broad_grouping,
      y = n_genes,
      fill = col_text
    )
  ) +
  geom_bar(
    stat = 'identity',
    colour = 'black'
  ) +
  scale_fill_manual(values = colPal(15)) +
  coord_flip() +
  shadowtext::geom_shadowtext(data = curated.manual.annotation[1,], mapping = aes(y = 17, label = col_text), colour = 'white', hjust = 0, fontface = 'bold') +
  shadowtext::geom_shadowtext(data = curated.manual.annotation[2,], mapping = aes(y = 4, label = col_text), colour = 'white', hjust = 0, fontface  = 'bold') +
  shadowtext::geom_shadowtext(data = curated.manual.annotation[3,], mapping = aes(y = 6, label = col_text), colour = 'white', hjust = 0, fontface  = 'bold') +
  shadowtext::geom_shadowtext(data = curated.manual.annotation[4,], mapping = aes(y = 5, label = col_text), colour = 'white', hjust = 0, fontface  = 'bold') +
  shadowtext::geom_shadowtext(data = curated.manual.annotation[5,], mapping = aes(y = 4, label = col_text), colour = 'white', hjust = 0, fontface  = 'bold') +
  geom_text(data = curated.manual.annotation[6,], mapping = aes(y = 20, label = col_text), colour = 'black', hjust = 0, fontface = 'bold') +
  geom_text(data = curated.manual.annotation[7,], mapping = aes(y = 17, label = col_text), colour = 'black', hjust = 0, fontface = 'bold') +
  geom_text(data = curated.manual.annotation[8,], mapping = aes(y = 15, label = col_text), colour = 'black', hjust = 0, fontface = 'bold') +
  geom_text(data = curated.manual.annotation[9,], mapping = aes(y = 15, label = col_text), colour = 'black', hjust = 0, fontface = 'bold') +
  geom_text(data = curated.manual.annotation[10,], mapping = aes(y = 12, label = col_text), colour = 'black', hjust = 0, fontface  = 'bold') +
  geom_text(data = curated.manual.annotation[11,], mapping = aes(y = 10, label = col_text), colour = 'black', hjust = 0, fontface  = 'bold') +
  geom_text(data = curated.manual.annotation[12,], mapping = aes(y = 10, label = col_text), colour = 'black', hjust = 0, fontface  = 'bold') +
  geom_text(data = curated.manual.annotation[13,], mapping = aes(y = 10, label = col_text), colour = 'black', hjust = 0, fontface  = 'bold') +
  geom_text(data = curated.manual.annotation[14,], mapping = aes(y = 8, label = col_text), colour = 'black', hjust = 0, fontface = 'bold') +
  geom_text(data = curated.manual.annotation[15,], mapping = aes(y = 8, label = col_text), colour = 'black', hjust = 0, fontface = 'bold') +
  scale_x_discrete(
    name = NULL,
    labels = rev(curated.manual.annotation$icon),
    limits = rev(levels(curated.manual.annotation$broad_grouping))
  ) +
  # ggtitle(label = 'Biological Groupings of Curated Genes<br>Under Selection in *A. laevis*') +
  theme_classic() +
  theme(
    # Axis
    axis.text.y = ggtext::element_markdown(color = 'black'),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    
    # Legend
    legend.position = 'none'
    
    # Plot title
    # plot.title = ggtext::element_markdown(face = 'bold', hjust = 0.5)
  )

barplot.psg.rate <- comb %>%
  arrange(rev(label)) %>%
  ggplot(aes(x = label, y = rate, fill = label)) +
  geom_bar(stat = 'identity', colour = 'black') +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(comb$label))) +
  shadowtext::geom_shadowtext(
    data = arrange(comb, rev(label))[1,],
    aes(label = rate, y = 5.5),
    fontface = 'bold',
    colour = 'white',
    size = 5
  ) +
  geom_text(
    data = arrange(comb, rev(label))[2,],
    aes(label = rate, y = 8),
    fontface = 'bold',
    colour = 'black',
    size = 5
  ) +
  shadowtext::geom_shadowtext(
    data = arrange(comb, rev(label))[3,],
    aes(label = rate, y = 10.5),
    fontface = 'bold',
    colour = 'white',
    size = 5
  ) +
  geom_text(
    data = arrange(comb, rev(label))[4,],
    aes(label = rate, y = 9),
    fontface = 'bold',
    colour = 'black',
    size = 5
  ) +
  labs(
    y = 'Number of PSGs'
    # title = 'Need to put\ntitle here'
  ) +
  scale_fill_manual(values = c('#2a9d8f', '#118ab2', '#ffb703', '#e76f51')) +
  theme_classic() +
  theme(
    # Elements that need to be blank
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    
    # Axis title/text
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    
    # Title
    # plot.title = element_text(face = 'bold', hjust = 0.5),
    
    # White space around plots
    plot.margin = margin(
      t = 1, 
      r = 0,
      b = 1, 
      l = 0, 
      unit = "cm"
    )
  )

treeplot <- ggtree(tree.data) +
  theme_tree2() +
  scale_x_continuous(
    labels = abs
  ) +
  xlab('\nMya')

treeplot <- revts(treeplot) +
  geom_label(aes(x = branch, label=label), fontface = 'bold.italic') +
  geom_tippoint(size = 3) +
  theme(
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14)
  )

# Dataset-4
full.plot <- (treeplot + barplot.psg.rate + plot_spacer() + barplot.annotation) + 
  plot_layout(widths = c(1, 0.3, 0.1, 1.5)) +
  plot_annotation(tag_levels = 'A')

ragg::agg_png(
  filename = here('manuscript', 'figures', 'figure-1.png'),
  height = 1000,
  width = 1500,
  units = 'px',
  res = 144
)
print(full.plot)
invisible(dev.off())
