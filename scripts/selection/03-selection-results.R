# ---------------------------------------------------------------------------- #
# Extract Selection Results
#
# - aBSREL
#   - Significant genes (p <= 0.001)
#   - Stratify into marine and terrestrial
#   - Within stratified groups, categorise significant genes by selection type
#     - Purifying (omega < 1)
#     - Positive-episodic (omega > 1)
#
# - RELAX
#   - Describe the type of natural selection occurring on the significant genes
#     - Intensification (K > 1)
#     - Relaxation (K < 1)
#
# - MEME (D9 only)
#   - Identify individual sites under selection
#     - Purifying (omega < 1)
#     - Positive-episodic (omega > 1)

# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(here)
})

# ---------------------------------------------------------------------------- #
# Import HyPhy data
hyphy_d4 <- read_rds(
  file = here('data', 'hyphy-dataframes', 'hyphy-four-sample.rds')
)

hyphy_d9 <- read_rds(
  file = here('data', 'hyphy-dataframes', 'hyphy-nine-sample.rds')
)

hyphy <- list(
  'Four' = hyphy_d4,
  'Nine' = hyphy_d9
)

rm(hyphy_d4)
rm(hyphy_d9)

gc()

genes_marine_D4 <- read_lines(
  here('data', 'selection-results', 'four', 'absrel-pval-0.001-marine.txt')
)
genes_marine_D9 <- read_lines(
  here('data', 'selection-results', 'nine', 'absrel-pval-0.001-marine.txt')
)

genes_ter_D4 <- read_lines(
  here('data', 'selection-results', 'four', 'absrel-pval-0.001-terrestrial.txt')
)
genes_ter_D9 <- read_lines(
  here('data', 'selection-results', 'nine', 'absrel-pval-0.001-terrestrial.txt')
)

annotation <- read_csv(
  file = here('data', 'utility-data', 'master-GO-annotation.csv'),
  col_names = TRUE, 
  col_types = cols()
)

# ---------------------------------------------------------------------------- #
# aBSREL

# Filter for p <= 0.001 and unnest rate classes
absrel_sig <- imap_dfr(hyphy, .id = 'Dataset', ~{
  
  if (.y == 'Four') {
    species <- 'aipysurusLaevis'
  } else {
    species <- 'Node_aipysurusLaevis'
  }
  
  # Gene length
  df_length <- .x[['absrel']][['Input']] %>%
    select(ID, `Number of Sites`)
  
  df <- .x[['absrel']][['Branch Attributes']] %>%
    filter(`Corrected P-value` <= 0.001) %>%
    left_join(df_length) %>%
    
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
        Species %in% species,
        'Marine',
        'Terrestrial'
      )
    )
}) %>%
  pivot_wider(names_from = value, values_from = `Rate Distribution`) %>%
  pivot_longer(
    18:23,
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

# Count histogram of proportion of sites belonging to each rate class
ragg::agg_png(
  filename = here('figures', 'absrel', 'histogram-proportion-sites-omega-rate-classes.png'),
  width = 2000,
  height = 1500,
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
    x = 'Proportion of Sites',
    y = 'Count',
    fill = expression("\u03C9 Rate Class")
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
    strip.text = element_text(size = 16, face = 'bold')
  ) +
  facet_grid(
    Condition ~ Dataset,
    scales = 'free_y'
  )
invisible(dev.off())

# Correlation of proportion of sites belonging to rate class against branch length
ragg::agg_png(
  filename = here('figures', 'absrel', 'point-proportion-of-sites-vs-full-model.png'),
  width = 2000,
  height = 1500,
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
    x = 'Full Adaptive Model (Branch lengths)',
    y = 'Proportion of Sites',
    colour = expression("\u03C9 Rate Class"),
    shape = expression("\u03C9 Rate Class")
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
    strip.text = element_text(size = 16, face = 'bold')
  ) +
  facet_grid(
    Condition ~ Dataset,
    scales = 'free_y'
  )
invisible(dev.off())

# Count histogram of omega estimates by rate class
ragg::agg_png(
  filename = here('figures', 'absrel', 'histogram-omega-values.png'),
  width = 2000,
  height = 1500,
  units = 'px',
  res = 144
)
absrel_sig %>%
  ggplot(aes(x = omega, fill = `Rate Class`)) +
  geom_histogram(
    position = 'identity',
    bins = 100
  ) +
  scale_fill_manual(
    values = c('#f4d75e', '#e9723d', '#0b7fab')
  ) +
  labs(
    x = expression("\u03C9"),
    y = 'Count',
    fill = expression("\u03C9 Rate Class")
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
    strip.text = element_text(size = 16, face = 'bold')
  ) +
  facet_grid(
    Condition ~ Dataset,
    scales = 'free'
  )
invisible(dev.off())

ragg::agg_png(
  filename = here('figures', 'absrel', 'point-length-omega-proportion.png'),
  width = 2000,
  height = 1500,
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
    x = 'Number of Sites',
    y = 'Proportion of Sites',
    colour = expression("\u03C9 Rate Class"),
    shape = expression("\u03C9 Rate Class")
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
    strip.text = element_text(size = 16, face = 'bold')
  ) +
  facet_grid(
    Condition ~ Dataset,
    scales = 'free_y'
  )
invisible(dev.off())

# ---------------------------------------------------------------------------- #
# RELAX
fs::dir_create(
  path = here('figures', 'relax')
)

relax_sig <- imap_dfr(hyphy, .id = 'Dataset', ~{
  
  if (.y == 'Four') {
    genes_marine <- genes_marine_D4
    genes_ter <- genes_ter_D4
  } else {
    genes_marine <- genes_marine_D9
    genes_ter <- genes_ter_D9
  }
  
  # Gene length
  df_length <- .x[['relax']][['Input']] %>%
    select(ID, `Number of Sites`)
  
  df <- .x[['relax']][['Test Results']] %>%
    # filter(`p-value` <= 0.001) %>%
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
          ID %in% genes_marine,
          'Marine',
          ifelse(
            ID %in% genes_ter,
            'Terrestrial',
            'Rest'
          )
        )
    )
}) %>%
  mutate(
    Group = factor(x = Group, levels = c('Marine', 'Terrestrial', 'Rest')),
    Significance = factor(x = Significance, levels = c('Significant', 'Insignificant'))
  )

# Sig/Insignificant - stratified by dataset and type of signal
ragg::agg_png(
  filename = here('figures', 'relax', 'bar-plot-summary-significant-signal.png'),
  width = 2000,
  height = 1500,
  units = 'px',
  res = 144
)
relax_sig %>%
  ggplot(aes(x = Condition, fill = Condition)) +
  geom_histogram(stat = 'count', colour = 'black') +
  scale_fill_manual(values = c('#F4CC70', '#DE7A22')) +
  theme_bw() +
  theme(
    # Legend
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = 'bottom',
    
    # Strip
    strip.text = element_text(size = 16, face = 'bold'),
    
    # Axis title/text
    axis.title = element_text(size = 16, face = 'bold'),
    axis.text = element_text(size = 14)
  ) +
  facet_grid(
    Significance ~ Dataset,
    scales = 'free_y'
  )
invisible(dev.off())

# Type of selection by group - marine/terrestrial genes
relax_sig %>%
  split(x = ., f = .$Dataset) %>%
  iwalk(~{
    path <- here('figures', 'relax', paste0(.y, '-bar-plot-significant-by-group.png'))
    ragg::agg_png(
      filename = path,
      width = 2000,
      height = 1500,
      units = 'px',
      res = 144
    )
    p <- .x %>%
      ggplot(aes(
        x = Condition,
        fill = Condition
      )) +
      geom_histogram(stat = 'count', colour = 'black') +
      scale_fill_manual(values = c('#F4CC70', '#DE7A22')) +
      theme_bw() +
      theme(
        # Legend
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = 'bottom',
        
        # Strip
        strip.text = element_text(size = 16, face = 'bold'),
        
        # Axis title/text
        axis.title = element_text(size = 16, face = 'bold'),
        axis.text = element_text(size = 14)
      ) +
      facet_grid(
        Significance ~ Group,
        scales = 'free_y'
      )
    print(p)
    invisible(dev.off())
  })


# ---------------------------------------------------------------------------- #
# MEME Plots - Unused
nine_sample_meme <- read_lines(
  paste0(
    'data/selection/nine-sample/genes-selection/nine-sample-absrel-pval-',
    pval,
    '-aipysurusLaevis.txt'
  )
)

sixteen_sample_meme <- read_lines(
  paste0('data/selection/sixteen-sample/genes-selection/sixteen-sample-absrel-pval-',
         pval,
         '-aipysurus.txt')
)

meme_list <- list(
  'nine-sample' = nine_sample_meme,
  'sixteen-sample' = sixteen_sample_meme
)

t <- imap(meme_list, ~{
  
  lst <- hyphy[[.y]]
  path <- file.path('data/figures', paste0(.y, '-pvalue-', pval, '-MEME.png'))
  a <- anno[ anno$ogid %in% .x, ] %>%
    distinct()
  
  df <- lst$meme$MLE %>%
    mutate(ID = sub('.*-', '', ID)) %>%
    filter(ID %in% .x) %>%
    left_join(a, by = c('ID' = 'ogid')) %>%
    mutate(ID = paste0(symbol, ' (', ID, ')')) %>%
    group_by(ID) %>%
    mutate(Position = 1:n()) %>%
    ungroup()
  
  p <- ggplot(data = df) +
    geom_point(data = filter(df, `p-value` > 0.1),
               aes(x = Position,
                   y = (1 - `p-value`),
                   colour = `p-value`),
               size = 3, colour = 'black') +
    geom_point(data = filter(df, `p-value` <= 0.05),
               aes(x = Position,
                   y = (1 - `p-value`),
                   colour = `p-value`),
               size = 3, colour = 'red') +
    geom_point(data = filter(df, `p-value` > 0.05, `p-value` <= 0.1),
               aes(x = Position,
                   y = (1 - `p-value`),
                   colour = `p-value`),
               size = 3, colour = RColorBrewer::brewer.pal(n = 3, name = 'Dark2')[1]) +
    labs(
      # y = expression(beta^'+'~'(Non-synonymous substitution rate for positive/neutral component at each site)'),
      y = 'Pvalue (1 - pvalue)',
      x = 'Position'
    ) +
    # guides(
    #   colour = guide_legend(title = expression('P-value'<='0.1'))
    # ) +
    theme_bw() +
    facet_wrap(. ~ ID, ncol = 3, scales = 'free_x')
  
  
  # p <- ggplot(data = df) +
  #   geom_point(data = filter(df, `p-value` > 0.1),
  #              aes(x = Position,
  #                  y = `beta+`,
  #                  colour = `p-value`),
  #              size = 3, colour = 'black') +
  #   geom_point(data = filter(df, `p-value` <= 0.05),
  #              aes(x = Position,
  #                  y = `beta+`,
  #                  colour = `p-value`),
  #              size = 3, colour = 'red') +
  #   geom_point(data = filter(df, `p-value` > 0.05, `p-value` <= 0.1),
  #              aes(x = Position,
  #                  y = `beta+`,
  #                  colour = `p-value`),
  #              size = 3, colour = RColorBrewer::brewer.pal(n = 3, name = 'Dark2')[1]) +
  #   labs(y = expression(beta^'+'~'(Non-synonymous substitution rate for positive/neutral component at each site)'),
  #        x = 'Position') +
  #   guides(colour = guide_legend(title = expression('P-value'<='0.1'))) +
  #   theme_bw() +
  #   facet_wrap(. ~ ID, ncol = 3, scales = 'free')
  ragg::agg_png(filename = path,
                width = 1200,
                height = 1200,
                units = 'px',
                res = 144)
  print(p)
  invisible(dev.off())
})

