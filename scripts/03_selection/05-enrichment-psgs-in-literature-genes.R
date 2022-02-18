# ------------------------------------------------------------------------------------------------ #
# Helpful links

# https://statsandr.com/blog/fisher-s-exact-test-in-r-independence-test-for-a-small-sample/
# https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/


# ------------------------------------------------------------------------------------------------ #
# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
})

# ------------------------------------------------------------------------------------------------ #
# Annotation object
annotation <- read_csv('data/utility-data/master-GO-annotation.csv.gz')

# ------------------------------------------------------------------------------------------------ #
# Import orthologs and get gene symbols
orthologs <- fs::dir_ls(
  path = 'data/oma-16-samples/complete-oma-groups/four-sample/nucleotide',
  all = TRUE
) %>%
  sub(
    '.*/(.*).fa', '\\1', .
  )

orthologs.genes <- annotation %>%
  filter(Group %in% orthologs) %>%
  select(Group, Gene) %>%
  distinct() %>%
  pull(Gene)

# ------------------------------------------------------------------------------------------------ #
# Import positively selected genes and get symbols
psgs <- psgs <- read_lines('data/selection-results/four/absrel-pval-0.001-marine.txt')
psgs.genes <- annotation %>%
  filter(Group %in% psgs) %>%
  pull(Gene) %>%
  sort() %>%
  unique()

# ------------------------------------------------------------------------------------------------ #
# Import literature genes (already symbols)
literature.genes <- read_csv(file = 'data/literature-marine-genes/literature-marine-genes.csv') %>%
  pull(symbol) %>%
  sort() %>%
  unique()

# ------------------------------------------------------------------------------------------------ #
# Is there significant overlap between our PSGs and marine literature

# Variables
group.psg = length(psgs.genes)
group.literature = length(literature.genes) - (length(intersect(literature.genes, orthologs.genes)))
overlap = length(intersect(literature.genes, psgs.genes))
total = (length(orthologs) + length(literature.genes)) - length(intersect(literature.genes, orthologs.genes))

# 2x2 contingency table
contigency.table <- matrix(
  data = c(
    overlap, 
    group.literature - overlap, 
    group.psg - overlap, total - group.literature - group.psg + overlap
  ),
  nrow = 2,
  dimnames = list(
    PSGs = c("Yes", "No"),
    Literature = c("Yes", "No")
  )
)

# Enrichment - pure p-value
result.phyper <- phyper(
  overlap - 1, 
  group.literature, 
  total - group.literature, 
  group.psg, 
  lower.tail = FALSE
)

# Fisher test - over-representation
result.fisher <- fisher.test(
  contigency.table,
  alternative='greater'
)
print(result.fisher)

# ------------------------------------------------------------------------------------------------ #
# Testing - only use the 157 marine genes that overlap our 3,594
#     - This means our gene universe stays as 3,594 as we're only considering marine genes
#       we could actually overlap with in the first place.
#     - The remaining 138 non-intersecting genes would be because they are not marine specific
test <- matrix(
  data = c(
    19, 
    157 - 19, 
    360, 
    3077
  ),
  nrow = 2,
  dimnames = list(
    PSGs = c("Yes", "No"),
    Literature = c("Yes", "No")
  )
)

fisher.test(
  test,
  alternative='greater'
)$p.value
