# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
})

# ---------------------------------------------------------------------------- #
# Data
gene.summary <- read_csv(
  file = here('manuscript', 'tables', 'aBSREL-p-value-summary.csv'),
  col_names = TRUE,
  col_types = cols()
)

# ---------------------------------------------------------------------------- #
# Plot
ragg::agg_png(
  filename = here('manuscript', 'figures', 'supplementary-aBSREL-pvalue-summary.png'),
  width = 1700,
  height = 1000,
  units = 'px',
  res = 144
)
gene.summary %>%
  pivot_longer(
    cols = 3:4,
    names_to = 'Dataset',
    values_to = 'Number of Orthologs',
  ) %>%
  mutate(
    `P-value` = factor(x = `P-value`, levels = c('0.05', '0.01', '0.001'))
  ) %>%
  filter(str_detect(string = Condition, pattern = 'unique')) %>%
  ggplot(
    aes(
      x = `P-value`,
      y = `Number of Orthologs`,
      fill = `P-value`
    )
  ) +
  geom_bar(
    alpha = 0.8,
    colour = 'black',
    position = 'dodge', stat = 'identity'
  ) +
  scale_fill_manual(values = wesanderson::wes_palette(name = 'Darjeeling1', n = 3)) + 
  # ggpomological::scale_fill_pomological() +
  theme_bw() +
  theme(
    # Legend
    legend.position = 'none',
    # legend.title = element_text(size = 16, face = 'bold'),
    # legend.text = element_text(size = 14),
    
    # Axes
    axis.title = element_text(size = 16, face = 'bold'),
    axis.text = element_text(size = 14),
    
    # Strip text
    strip.text.x = element_text(size = 14, face = 'bold'),
    strip.text.y = element_text(size = 14, face = 'bold')
  ) +
  facet_grid(
    Dataset ~ Condition
    # scales = 'free_y'
  )
invisible(dev.off())
 
