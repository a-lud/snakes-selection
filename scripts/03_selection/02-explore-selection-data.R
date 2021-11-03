# ---------------------------------------------------------------------------- #
# Explore the selection results

# ---------------------------------------------------------------------------- #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
})

# ---------------------------------------------------------------------------- #
# Import HyPhy data
hyphy_d4 <- read_rds(
  file = 'data/hyphy-dataframes/hyphy-four-sample.rds'
)

hyphy_d9 <- read_rds(
  file = 'data/hyphy-dataframes/hyphy-nine-sample.rds'
)

# List object to iterate over
hyphy <- list(
  'four' = hyphy_d4,
  'nine' = hyphy_d9
)

rm(hyphy_d4)
rm(hyphy_d9)

gc()

# ---------------------------------------------------------------------------- #
# Corrected P-value
pval <- list(
  '05' = 0.05, 
  '01' = 0.01, 
  '001' = 0.001
)

# ---------------------------------------------------------------------------- #
# Annotation object
anno <- read_rds('data/utility-data/master-annotation.rds') %>%
  rename(symbol = Gene, ogid = Group) %>%
  select(ogid, symbol)

# ---------------------------------------------------------------------------- #
# Positive test results
iwalk(hyphy, ~{
  # Output directory
  fs::dir_create(
    path = 'figures/absrel',
    recurse = TRUE
  )
  
  # Output file
  path <- file.path(
    'figures/absrel', 
    paste0(.y, '-positive-test-results.png')
  )
  
  # Plot
  p <- .x$absrel$`Test results` %>%
    ggplot(aes(x = `Positive test results`)) +
    geom_bar(
      aes(y = ..prop..),
      colour = 'black', 
      fill = 'grey',
      alpha = 0.6
    ) +
    scale_y_continuous(
      name = 'Proportion of Genes', 
      breaks = seq(0, 0.7, 0.1)
    ) +
    scale_x_continuous(name = 'Number of Branches With a Signal of Selection') +
    theme_bw() +
    theme(
      axis.title = element_text(size = 16, face = 'bold'),
      axis.text = element_text(size = 14)
    )
  ragg::agg_png(filename = path, 
                width = 1000,
                height = 1000, 
                units = 'px', 
                res = 144)
  print(p)
  invisible(dev.off())
})

# ---------------------------------------------------------------------------- #
# Number of Rate Classes
iwalk(hyphy, ~{
  
  # Output file
  path <- file.path(
    'figures/absrel', 
    paste0(.y, '-rate-classes.png')
  )
  
  # Get 'foreground' species
  test_branches <- .x$absrel$Tested %>% 
    filter(Status == 'test') %>% 
    extract2('Node') %>% 
    unique()
  
  p <- .x$absrel$`Branch Attributes` %>%
    filter(Species %in% test_branches) %>%
    mutate(
      Species = case_when(
        Species == 'aipysurusLaevis' ~ 'Aipysurus laevis',
        Species == 'najaNaja' ~ 'Naja naja',
        Species == 'notechisScutatus' ~ 'Notechis scutatus',
        Species == 'pseudonajaTextilis' ~ 'Pseudonaja textilis',
        Species == 'Node_aipysurusLaevis' ~ 'Node Aipysurus laevis',
        Species == 'Node_notechisScutatus' ~ 'Node Notechis scutatus',
        Species == 'Node_pseudonajaTextilis' ~ 'Node Pseudonaja textilis',
      )
    ) %>%
    ggplot(aes(x = `Rate Classes`)) +
    geom_bar(
      aes(y = ..prop..), 
      colour = 'black', 
      fill = 'grey',
      alpha = 0.6
    ) +
    facet_wrap(.~Species) +
    scale_y_continuous(
      name = 'Proportion of Orthologues', 
      breaks = seq(0,1,0.2)
    ) +
    scale_x_continuous(name = 'Rate Class Categories') +
    theme_bw() +
    theme(
      axis.title = element_text(size = 16, face = 'bold'),
      axis.text = element_text(size = 14),
      strip.text.x = element_text(size = 14, face = 'bold.italic')
    )
  ragg::agg_png(filename = path, 
                width = 1000,
                height = 1000, 
                units = 'px', 
                res = 144)
  print(p)
  invisible(dev.off())
})

 # ---------------------------------------------------------------------------- #
# Selection Summary Table: Number of genes passing alpha <= pval
imap(hyphy, ~{
  # Create output directory
  fs::dir_create(
    path = paste0('data/selection-results/', .y),
    recurse = TRUE
  )
  
  # Iterate over the p-values
  map(pval, function(p) {
    df <- .x[['absrel']][['Branch Attributes']] %>%
      filter(`Corrected P-value` <= p)
    
    # Specify marine node based on dataset
    if (.y == 'four') {
      mar <- 'aipysurusLaevis'
    } else {
      mar <- 'Node_aipysurusLaevis'
    }
    
    # Build summary table
    marine <- filter(df, Species == mar) %>% 
      distinct(ID) %>% 
      extract2('ID') %>%
      gsub('.*-', '', .) %>%
      unique()
    
    terrestrial <- filter(df, Species != mar) %>%
      distinct(ID) %>%
      extract2('ID') %>%
      gsub('.*-', '', .) %>%
      unique()
    
    write(x = setdiff(marine, terrestrial),
          file = file.path(
            'data/selection-results', 
            .y, 
            paste0('absrel-pval-', p, '-marine.txt'))
    )
    write(x = setdiff(terrestrial, marine),
          file = file.path(
            'data/selection-results', 
            .y,
            paste0('absrel-pval-', p, '-terrestrial.txt'))
    )
    
    df_return <- tibble::tribble(
      ~Condition, ~`Number of genes`,
      "aBSREL Aipysurus laevis", length(marine),
      "aBSREL Aipysurus laevis unique", length(setdiff(marine, terrestrial)),
      "aBSREL terrestrial", length(terrestrial),
      "aBSREL terrestrial unique", length(setdiff(terrestrial, marine))
    ) %>% mutate(
      dataset = .y,
      threshold = p
    )
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  pivot_wider(names_from = dataset, values_from = `Number of genes`) %>%
  mutate(
    Condition = factor(
      Condition,
      levels = c('aBSREL Aipysurus laevis',
                 'aBSREL Aipysurus laevis unique',
                 'aBSREL terrestrial',
                 'aBSREL terrestrial unique')
    ),
    threshold = factor(
      threshold,
      levels = c(
        0.05,
        0.01,
        0.001
      )
    )
  ) %>%
  arrange(threshold, Condition) %>%
  write_csv(
    x = ., 
    file = 'data/selection-results/absrel-summary.csv'
  )

# ---------------------------------------------------------------------------- #
# Venn Diagram: Selection Overlap
ogids <- imap(hyphy, ~{
  # Node to filter on to get marine/terrestrial
  if(.y == 'four') {
    string <- 'aipysurusLaevis'
  } else {
    string <- 'Node_aipysurusLaevis'
  }
  
  # Plot VennDiagram and return vector of OGIDs
  ogids <- map(pval, function(p) {
    # File path for output
    path <- file.path(
      'figures',
      paste0(
        .y, 
        '-venn-absrel-relax-',
        p,
        '.png'
      )
    )
    
    # Get sites passing thresholds
    m <- .x$absrel$`Branch Attributes` %>%
      filter(`Corrected P-value` <= p, Species == string) %>%
      mutate(ID = sub('.*-', '', ID)) %>%
      extract2('ID') %>%
      unique()
    
    t <- .x$absrel$`Branch Attributes` %>%
      filter(`Corrected P-value` <= p, Species != string) %>%
      mutate(ID = sub('.*-', '', ID)) %>%
      extract2('ID') %>%
      unique()
    
    r <- .x$relax$`Test Results` %>%
      filter(`p-value` <= p) %>%
      mutate(ID = sub('.*-', '', ID)) %>%
      extract2('ID') %>%
      unique()
    
    silence <- VennDiagram::venn.diagram(
      x = list(r, m, t),
      filename = path,
      height = 1500,
      width = 1500,
      resolution = 400,
      imagetype = 'png',
      output = TRUE,
      
      # Categories
      category.names = c('RELAX', 'Marine-aBSREL', 'Terrestrial-aBSREL'),
      
      # Circles
      lwd = 2, # Line width of circles
      col = RColorBrewer::brewer.pal(n = 3, name = 'Set1'), # Colour of each circles circumference
      fill = RColorBrewer::brewer.pal(n = 3, name = 'Set1'), # Colour of each circle
      alpha = 0.5, # Alpha of colours
      
      # Numbers
      cex = 0.6, # Number size
      fontface = 'bold',
      fontfamily = 'sans',
      
      # Category names
      cat.pos = c(-27, 27, 135), # Category position
      cat.dist = c(0.055, 0.055, 0.085), # Category distance
      cat.cex = 0.6, # Category size
      cat.fontface = 'bold',
      cat.fontfamily = 'sans',
      cat.default.pos = 'outer',
      rotation = 1
    )
    
    list(
      'absrel-marine' = m, 
      'absrel-terrestrial' = t, 
      'relax' = r
    )
  })
})

# ---------------------------------------------------------------------------- #
# Overlap between four and nine sample datasets
df_intersect <- map2(.x = names(ogids$four), .y = names(ogids$nine), ~{
  
  four <- ogids$four[[.x]]
  nine <- ogids$nine[[.y]]
  
  # D4
  tibble(
    'Model' = c('aBSREL-Marine', 'aBSREL-Terrestrial', 'RELAX'),
    'D4' = c(
      length(four[['absrel-marine']]), 
      length(four[['absrel-terrestrial']]), 
      length(four[['relax']])
    ),
    'D9' = c(
      length(nine[['absrel-marine']]), 
      length(nine[['absrel-terrestrial']]), 
      length(nine[['relax']])
    ),
    'Intersect' = c(
      length(intersect(four[['absrel-marine']], nine[['absrel-marine']])),
      length(intersect(four[['absrel-terrestrial']], nine[['absrel-terrestrial']])),
      length(intersect(four[['relax']], nine[['relax']]))
    ))
}) %>%
  set_names(value = pval) %>%
  bind_rows(.id = 'threshold') %>%
  mutate(threshold = as.double(threshold))

write_csv(
  x = df_intersect,
  file = 'data/selection-results/intersect-D4-D9-results.csv',
  col_names = TRUE
)

# ---------------------------------------------------------------------------- #
# Upset plots - TODO: Italicise species
combination_matrix <- imap(hyphy, ~{
  map(pval, function(p) {
    path <- file.path(
      'figures',
      paste0(.y, '-upset-', p, '.png')
    )
    
    # Set labels
    if(.y == 'four'){
      string1 <- 'aipysurusLaevis'
      string2 <- 'notechisScutatus'
      string3 <- 'pseudonajaTextilis'
      string4 <- 'najaNaja'
    } else {
      string1 <- 'Node_aipysurusLaevis'
      string2 <- 'Node_notechisScutatus'
      string3 <- 'Node_pseudonajaTextilis'
      string4 <- 'najaNaja'
    }
    
    # Extract OGIDs
    s1 <- .x$absrel$`Branch Attributes` %>%
      filter(Species == string1,
             `Corrected P-value` <= p) %>%
      mutate(ID = sub('.*-', '', ID)) %>%
      extract2('ID')
    s2 <- .x$absrel$`Branch Attributes` %>%
      filter(Species == string2,
             `Corrected P-value` <= p) %>%
      mutate(ID = sub('.*-', '', ID)) %>%
      extract2('ID')
    s3 <- .x$absrel$`Branch Attributes` %>%
      filter(Species == string3,
             `Corrected P-value` <= p) %>%
      mutate(ID = sub('.*-', '', ID)) %>%
      extract2('ID')
    s4 <- .x$absrel$`Branch Attributes` %>%
      filter(Species == string4,
             `Corrected P-value` <= p) %>%
      mutate(ID = sub('.*-', '', ID)) %>%
      extract2('ID')
    
    # Named list for ComplexHeatmap's UpSet
    l <- list(
      'Aipysurus laevis' = s1,
      'Notechis scutatus' = s2,
      'Pseudonaja textilis' = s3,
      'Naja naja' = s4
    )
    
    # Combinations matrix
    mat <- ComplexHeatmap::make_comb_mat(l)
    
    # Upset plot
    upSet <- ComplexHeatmap::UpSet(
      m = mat,
      pt_size = unit(5, "mm"), 
      lwd = 3,
      set_order = names(sort(ComplexHeatmap::set_size(mat))),
      right_annotation = NULL,
      left_annotation = ComplexHeatmap::rowAnnotation(
        "Set size" = ComplexHeatmap::anno_barplot(
          axis_param = list(
            gp = grid::gpar(fontsize=12)
          ),
          ComplexHeatmap::set_size(mat), 
          border = FALSE,
          gp = grid::gpar(fill = "black"), 
          width = unit(3, "cm")
        )
      ),
      comb_col = c("#c03728", "#919c4c", "#f5c04a", "#828585")[ComplexHeatmap::comb_degree(mat)],
      # column_title = "UpSet: Common Orthologs Between\naBSREL-Marine and aBSREL-Terrestrial",
      row_names_gp = grid::gpar(fontface = 'italic')
    )
    
    ragg::agg_png(filename = path,
                  width = 1000,
                  height = 1000,
                  units = 'px',
                  res = 144)
    print(upSet)
    invisible(dev.off())

    # Return combination matrix
    return(mat)
  })
})

# ---------------------------------------------------------------------------- #
# Combination matrices
iwalk(combination_matrix, ~{
  imap(.x, function(mat, threshold) {
    tibble(
      id = names(ComplexHeatmap::comb_size(mat)),
      intersect = ComplexHeatmap::comb_size(mat)
    ) %>%
      separate(
        col = id,
        into = c('A. laevis', 'N. scutatus', 'P. textilis', 'N. naja'),
        sep = "(?<=.)"
      ) %>%
      arrange(-intersect) %>%
      select(Intersect = intersect, everything())
  })
}) %>%
  write_rds(
    x = ., 
    file = 'data/selection-results/combination-matrices.rds'
  )
