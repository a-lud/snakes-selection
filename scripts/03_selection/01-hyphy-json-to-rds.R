# Convert HyPhy JSON files to Dataframes

# ---------------------------------------------------------------------------- #
# Scripts
source('scripts/util/updatedHyphyParser.R')

# ---------------------------------------------------------------------------- #
# Pre-run results
fs::dir_create(
  path = 'data/hyphy-dataframes', 
  recurse = TRUE
)

# ---------------------------------------------------------------------------- #
# D4 - Four snake dataset
results <- parseHyphy(
  path = 'data/hyphy/four-sample',
  analysis = c('absrel', 'busted', 'relax')
)

readr::write_rds(
  x = results,
  file = 'data/hyphy-dataframes/hyphy-four-sample.rds',
  compress = 'gz'
)

# ---------------------------------------------------------------------------- #
# Nine snake dataset
results <- parseHyphy(
  path = 'data/hyphy/nine-sample',
  analysis = c('busted', 'absrel', 'meme', 'relax')
)

readr::write_rds(
  x = results,
  file = 'data/hyphy-dataframes/hyphy-nine-sample.rds',
  compress = 'gz'
)
