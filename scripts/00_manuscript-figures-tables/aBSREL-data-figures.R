# ---------------------------------------------------------------------------- #
# aBSREL data cleaning + figures
#
# Generate summary table and figures relating to the aBSREL analysis. Figures
# include:
#   - Bar plot of orthologues with each omega rate class
#   - Histogram showing the proportion of sites belonging to each rate class
#   - Correlation of branch length against proportion of sites belonging to each
#     rate class
#   - Correlatino of gene length against proportion of sites belonging to each
#     rate class

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

fs::dir_create(
  path = here('manuscript', 'figures')
)

# ---------------------------------------------------------------------------- #
# Import data
hyphy <- read_rds(
  file = here('data', 'hyphy-dataframes', 'hyphy-four-sample.rds')
)

# Gene lengths
df_length <- hyphy[['absrel']][['Input']] %>%
  select(ID, `Number of Sites`)

# ---------------------------------------------------------------------------- #
# aBSREL
absrel_sig <- hyphy[['absrel']][['Branch Attributes']] %>%
  left_join(df_length) %>%
  mutate(ID = sub('.*-', '', ID)) %>%
  filter(`Corrected P-value` <= 0.001) %>%
  
  # Unnest rate classes and label accordingly (omega1, prop1, omega2, prop2, ...)
  unnest(cols = `Rate Distribution`) %>%
  group_by(ID, Species) %>%
  mutate(
    value = paste0(
      c('omega', 'proportion'),
      rep(seq(1:(n()/2)), each = 2)
    )
  ) %>%
  ungroup() %>%
  mutate(
    Condition = ifelse(
      Species == 'aipysurusLaevis',
      'Aipysurus laevis',
      'Terrestrial'
    )
  ) %>%
  pivot_wider(names_from = value, values_from = `Rate Distribution`) %>%
  pivot_longer(
    17:22,
    names_to = c('.value', 'Rate Class'),
    names_pattern = "(.*)(1|2|3)"
  ) %>%
  # Remove rows where omega and proportion are both NA (no need to keep them)
  filter(
    across(
      .cols = c('omega', 'proportion'),
      ~!is.na(.x)
    )
  )

# Write data to file
write_csv(
  x = absrel_sig,
  file = here('manuscript', 'tables', 'aBSREL-significant-table.csv'), 
  col_names = TRUE
)

# ---------------------------------------------------------------------------- #
# Figures

# Number of Rate Classes ---
test_branches <- hyphy$absrel$Tested %>% 
  filter(Status == 'test') %>% 
  extract2('Node') %>% 
  unique()

ragg::agg_png(
  filename = here('manuscript', 'figures',  'absrel-rate-classes.png'), 
  width = 1000,
  height = 1000, 
  units = 'px', 
  res = 144
)
hyphy$absrel$`Branch Attributes` %>%
  filter(Species %in% test_branches) %>%
  mutate(
    Species = case_when(
      Species == 'aipysurusLaevis' ~ 'Aipysurus laevis',
      Species == 'najaNaja' ~ 'Naja naja',
      Species == 'notechisScutatus' ~ 'Notechis scutatus',
      Species == 'pseudonajaTextilis' ~ 'Pseudonaja textilis'
    ),
    Species = factor(x = Species, levels = c('Aipysurus laevis', 'Notechis scutatus', 'Pseudonaja textilis', 'Naja naja'))
  ) %>%
  ggplot(aes(x = `Rate Classes`)) +
  geom_bar(
    aes(y = ..prop..), 
    colour = 'black', 
    fill = 'grey',
    alpha = 0.6
  ) +
  facet_wrap( . ~ Species) +
  scale_x_continuous(name = expression(bold("Omega (\u03C9) Rate Classes")))+
  scale_y_continuous(
    name = 'Proportion of Orthologues', 
    breaks = seq(0,1,0.2)
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16, face = 'bold'),
    axis.text = element_text(size = 14),
    strip.text.x = element_text(size = 14, face = 'bold.italic')
  )
invisible(dev.off())

# Proportion of sites along each gene belonging to each omega-rate-class ---
ragg::agg_png(
  filename = here('manuscript', 'figures', 'absrel-histogram-proportion-sites-omega-rate-classes.png'),
  width = 1500,
  height = 1000,
  units = 'px',
  res = 144
)
absrel_sig %>% 
  ggplot(aes(x = proportion, fill = `Rate Class`)) +
  geom_histogram(
    position = 'identity',
    bins = 100,
    # colour = 'black',
    alpha = 0.7
  ) +
  scale_fill_manual(
    values = c('#f4d75e', '#e9723d', '#0b7fab')
  ) +
  labs(
    x = '\nProportion of Sites Belonging\nto Each Rate Class',
    y = 'Count',
    fill = expression("Omega (\u03C9)\nRate Class")
  ) +
  theme_bw() +
  theme(
    # Legend
    legend.title = element_text(size = 16, face = 'bold'),
    legend.text = element_text(size = 14),
    
    # Axes
    axis.title = element_text(size = 16, face = 'bold'),
    axis.text = element_text(size = 14),
    
    # Strip text
    strip.text = element_text(size = 16, face = 'bold.italic')
  ) +
  facet_wrap(
    Condition ~ ., 
    scales = 'free_y',
    nrow = 2
  )
invisible(dev.off())

# Correlate proportion of sites belonging to rate class against branch length ---
ragg::agg_png(
  filename = here('manuscript', 'figures', 'absrel-correlate-branchLen-propSites.png'),
  width = 1500,
  height = 1000,
  units = 'px',
  res = 144
)
absrel_sig %>% 
  ggplot(aes(
    y = proportion,
    x = `Full adaptive model`, 
    colour = `Rate Class`,
    shape = `Rate Class`
  )) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(
    values = c('#f4d75e', '#e9723d', '#0b7fab')
  ) +
  labs(
    x = '\nFull Adaptive Model (Branch lengths)',
    y = 'Proportion of Sites Belonging\nto Each Rate Class\n',
    colour = expression("Omega (\u03C9)\nRate Class"),
    shape = expression("Omega (\u03C9)\nRate Class")
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 5))
  ) +
  theme_bw() +
  theme(
    # Legend
    legend.title = element_text(size = 16, face = 'bold'),
    legend.text = element_text(size = 14),
    
    # Axes
    axis.title = element_text(size = 16, face = 'bold'),
    axis.text = element_text(size = 14),
    
    # Strip text
    strip.text = element_text(size = 16, face = 'bold.italic')
  ) +
  facet_wrap(
    Condition ~ .,
    scales = 'free_y',
    nrow = 2
  )
invisible(dev.off())

# Length of gene vs proportion ---
ragg::agg_png(
  filename = here('manuscript', 'figures', 'absrel-correlate-length-omega-proportion.png'),
  width = 1500,
  height = 1000,
  units = 'px',
  res = 144
)
absrel_sig %>%
  ggplot(aes(
    x = `Number of Sites`,
    y = proportion,
    colour = `Rate Class`,
    shape = `Rate Class`
  )) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(
    values = c('#f4d75e', '#e9723d', '#0b7fab')
  ) +
  labs(
    x = '\nGene Length (bp)',
    y = 'Proportion of Sites Belonging\nto Each Rate Class\n',
    colour = expression("Omega (\u03C9)\nRate Class"),
    shape = expression("Omega (\u03C9)\nRate Class")
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 5))
  ) +
  theme_bw() +
  theme(
    # Legend
    legend.title = element_text(size = 16, face = 'bold'),
    legend.text = element_text(size = 14),
    
    # Axes
    axis.title = element_text(size = 16, face = 'bold'),
    axis.text = element_text(size = 14),
    
    # Strip text
    strip.text = element_text(size = 16, face = 'bold.italic')
  ) +
  facet_wrap(
    Condition ~ .,
    scales = 'free_y',
    nrow = 2
  )
invisible(dev.off())

