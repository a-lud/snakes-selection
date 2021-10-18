# ---------------------------------------------------------------------------- #
# Parse HyPhy JSON files

# These are quite slow for many files, so there are efficiency gains to be made.
# ---------------------------------------------------------------------------- #

# Main function
parseHyphy <- function(path, 
                       analysis = NULL) {
  
  # List all Json files
  files <- fs::dir_info(
    path = path, 
    recurse = TRUE, 
    type = 'file', 
    glob = '*.json'
  )
  
  files <- dplyr::filter(
    files,
    size != 0
  )
  
  files <- as.character(files[['path']])
  
  # Subset for analyeses of interest
  if (!is.null(analysis)) {
    pattern <- paste(analysis, collapse = '|')
  } else {
    analysis <- c('absrel', 'busted', 'meme', 'relax')
    pattern <- paste(analysis, collapse = '|')
  }
  
  files <- files[stringr::str_detect(
    string = files,
    pattern = pattern
  )]
  
  # Set vector names to match model-OGID
  names(files) <- sub('^(.*)-(OG\\d+)(?:-|_).*', '\\1-\\2', basename(files))
  
  # Split into list based on analysis
  files <- purrr::map(analysis, ~{
    files[stringr::str_detect(
      string = files,
      pattern = .x
    )]
  })
  
  names(files) <- analysis
  
  # Iterate over analyses and subset json files
  out <- purrr::imap(files, ~{
    lst_json <- purrr::map(.x, jsonlite::read_json)
    
    # Common fields
    input <- getInputs(lst_json)
    tested <- getTested(lst_json)
    dataPartition <- getDataPartitions(lst_json)
    
    if (.y == 'absrel') {
      print('Import: aBSREL')
      list(
        'Input' = input,
        'Tested' = tested,
        'Test results' = getTestResultsAbsrel(lst_json),
        'Data Partition' = dataPartition,
        'Branch Attributes' = getBranchAttributesAbsrel(lst_json),
        'Fits' = getFitsAbsrel(lst_json)
      )
    } else if (.y == 'busted') {
      print('Import: BUSTED')
      list(
        'Input' = input,
        'Tested' = tested,
        'Evidence Ratios' = getEvidenceRatiosBusted(lst_json),
        'Site Log-likelihoods' = getSiteLogLikelihoodBusted(lst_json),
        'Branch Attributes' = getBranchAttributesBusted(lst_json),
        'Data Partitions' = dataPartition,
        'Fits' = getFitsBusted(lst_json),
        'Test Results' = getTestResultsBusted(lst_json)
      )
    } else if (.y == 'meme') {
      print('Import: MEME')
      list(
        'MLE' = getMleMeme(lst_json),
        'Branch Attributes' = getBranchAttributesMeme(lst_json),
        'Data Partitions' = dataPartition,
        'Fits' = getFitsMeme(lst_json),
        'Input' = input,
        'Tested' = tested
      )
    } else if (.y == 'relax') {
      print('Import: RELAX')
      list(
        'Input' = input,
        'Tested' = tested,
        'Branch Atrributes' = getBranchAttributesRelax(lst_json),
        'Data partitions' = dataPartition,
        'Fits' = getFitsRelax(lst_json),
        'Test Results' = getTestResultsRelax(lst_json)
      )
    }
  })
  
  # Name output list + return
  return(out)
}

# ---------------------------------------------------------------------------- #
# Shared fields

# Get 'input' key-values
getInputs <- function(l) {
  
  out <- purrr::map(l, function(json) {
    input <- json[['input']]
    
    tibble::tibble(
      'File Name' = basename(input[['file name']]),
      'Number of Sequences' = as.integer(input[['number of sequences']]),
      'Number of Sites' = as.integer(input[['number of sites']]),
      'Partition Count' = as.integer(input[['partition count']]),
      'Trees' = list(unlist(input[['trees']]))
    )
  })
  
  return(dplyr::bind_rows(out, .id = 'ID'))
}

# Get 'tested' key-values
getTested <- function(l) {
  out <- purrr::map(l, function(json) {
    partitions <- json[['tested']]
    
    df <- purrr::map(partitions, function(pt) {
      tibble::tibble(
        'Node' = names(pt),
        'Status' = unlist(pt)
      )
    })
    dplyr::bind_rows(df, .id = 'Partition')
  })
  return(dplyr::bind_rows(out, .id = 'ID'))
}

# Get 'data partitions' key-values
getDataPartitions <- function(l) {
  out <- purrr::map(l, function(json) {
    partitions <- json[['data partitions']]
    
    df <- purrr::map(partitions, function(pt) {
      tibble::tibble(
        'Name' = pt[['name']],
        'Coverage' = list(unlist(pt[['coverage']]))
      )
    })
    dplyr::bind_rows(df, .id = 'Partition')
  })
  return(dplyr::bind_rows(out, .id = 'ID'))
}

# ---------------------------------------------------------------------------- #
# aBSREL specific

getTestResultsAbsrel <- function(l) {
  out <- purrr::map(l, function(json) {
    tibble::tibble(
      'P-value threshold' = json[['test results']][['P-value threshold']],
      'Positive test results' = json[['test results']][['positive test results']],
      'Tested' = json[['test results']][['tested']]
    )
  })
  
  return(dplyr::bind_rows(out, .id = 'ID'))
}

# Branch Attributes:
getBranchAttributesAbsrel <- function(l) {
  out <- purrr::map(l, function(json) {
    
    # Drop 'attributes' key from partitions (not needed)
    partitions <- json[['branch attributes']]
    partitions <- partitions[-length(partitions)]
    
    df_partition <- purrr::map(partitions, function(pt) {
      
      df_sp <- purrr::map(pt, function(sp) {
        tibble::tibble(
          # 'Species' = sp,
          'Baseline MG94xREV' = sp[['Baseline MG94xREV']],
          'Baseline MG94xREV omega ratio' = sp[['Baseline MG94xREV omega ratio']],
          'Corrected P-value' = ifelse(is.null(sp[['Corrected P-value']]), 
                                       NA_integer_, 
                                       sp[['Corrected P-value']]),
          'Full adaptive model' = sp[['Full adaptive model']],
          'Full adaptive model (non-synonymous subs/site)' = sp[['Full adaptive model (non-synonymous subs/site)']],
          'Full adaptive model (synonymous subs/site)' = sp[['Full adaptive model (synonymous subs/site)']],
          'LRT' = ifelse(is.null(sp[['LRT']]), NA_integer_, sp[['LRT']]),
          'Nucleotide GTR' = sp[['Nucleotide GTR']],
          'Rate Distribution' = list(unlist(sp[['Rate Distributions']])),
          'Rate Classes' = as.integer(sp[['Rate classes']]),
          'Uncorrected P-value' = ifelse(is.null(sp[['Uncorrected P-value']]), 
                                         NA_integer_, 
                                         sp[['Uncorrected P-value']]),
          'Original Name' = sp['Original Name']
        )
      })
      
      # Tibble of all species in analysis
      dplyr::bind_rows(df_sp, .id = 'Species')
    })
    
    # Tibble of all partitions in an analysis
    dplyr::bind_rows(df_partition, .id = 'Partition')
  })
  
  # Tibble across all sequences analysed
  return(dplyr::bind_rows(out, .id = 'ID'))
}

# Fits: MG94xREV, Full adaptive model, Nucleotide GTR
getFitsAbsrel <- function(l) {
  out <- purrr::map(l, function(json) {
    
    # Store model names and models into objects
    nms <- names(json[['fits']])
    mdl <- json[['fits']]
    
    # Iterate over each model and extract relevant information
    df_models <- purrr::map(nms, function(m) {
      if (m == 'Baseline MG94xREV') {
        # Get Alpha-sorted codons (excluding STOP codons)
        codons <- sort(names(Biostrings::getGeneticCode()))
        codons <- codons[ ! codons %in% c('TAG', 'TAA', 'TGA') ]
        
        # Get MG94xREV attributes
        aicc <- mdl[[m]][['AIC-c']]
        eqFreq <- magrittr::set_names(x = unlist(mdl[[m]][['Equilibrium frequencies']]),
                                      value = codons)
        loglike <- mdl[[m]][['Log Likelihood']]
        rateDist <- unlist(mdl[[m]][['Rate Distributions']][['Per-branch omega']])
        estPar <- mdl[[m]][['estimated parameters']]
        
        tibble::tibble(
          'AIC-c' = aicc,
          'Equilibrium Frequency' = list(eqFreq),
          'Log-likelihood' = loglike,
          'Rate Distribution' = list(rateDist),
          'Estimated parameters' = estPar
        )
      } else if (m == 'Full adaptive model') {
        aicc <- mdl[[m]][['AIC-c']]
        loglike <- mdl[[m]][['Log Likelihood']]
        estPar <- mdl[[m]][['estimated parameters']]
        
        tibble::tibble(
          'AIC-c' = aicc,
          'Log-likelihood' = loglike,
          'Estimated parameters' = estPar
        )
      } else {
        aicc <- mdl[[m]][['AIC-c']]
        eqFreq <-magrittr::set_names(x = unlist(mdl[[m]][['Equilibrium frequencies']]),
                                     value = c('A', 'C', 'G', 'T'))
        loglike <- mdl[[m]][['Log Likelihood']]
        rateDist <- unlist(mdl[[m]][['Rate Distributions']])
        estPar <- mdl[[m]][['estimated parameters']]
        
        tibble::tibble(
          'AIC-c' = aicc,
          'Equilibrium Frequency' = list(eqFreq),
          'Log-likelihood' = loglike,
          'Rate Distribution' = list(rateDist),
          'Estimated parameters' = estPar
        )
      }
    })
    
    # Bind results of each model per sequence
    names(df_models) <- nms
    df_models <- dplyr::bind_rows(df_models, .id = 'Model')
    df_models <- dplyr::arrange(.data = df_models,
                                dplyr::desc(`Log-likelihood`))
  })
  
  # Return dataframe of all sequences run
  return(dplyr::bind_rows(out, .id = 'ID'))
}

# ---------------------------------------------------------------------------- #
# BUSTED specific

# Get test-results
getTestResultsBusted <- function(l) {
  out <- purrr::map(l, function(json) {
    tibble::tibble(
      'LRT' = json[['test results']][['LRT']],
      'p-value' = json[['test results']][['p-value']],
    )
  })
  
  return(dplyr::bind_rows(out, .id = 'ID'))
}

# Get BUSTED fits data: 17.979
getFitsBusted <- function(l) {
  out <- purrr::map(l, function(json) {
    
    # Store model names and models into objects
    nms <- names(json[['fits']])
    mdl <- json[['fits']]
    
    # Iterate over each model and extract relevant information
    df_models <- purrr::map(nms, function(m) {
      if (m == 'Constrained model') {
        
        # Get MG94xREV attributes
        aicc <- mdl[[m]][['AIC-c']]
        loglike <- mdl[[m]][['Log Likelihood']]
        
        rateDist <- unlist(mdl[[m]][['Rate Distributions']])
        rateDist <- tibble::tibble(Names = names(rateDist), 
                                   Values = rateDist)
        rateDist <- tidyr::separate(data = rateDist,
                                    col = Names, 
                                    into = c('Branch set',
                                             'Omega rate classes',
                                             'propRate'), 
                                    sep = '\\.')
        rateDist <- dplyr::filter(.data = rateDist,
                                  `Branch set` == 'Test')
        rateDist <- tidyr::pivot_wider(data = rateDist,
                                       names_from = propRate,
                                       values_from = Values)
        rateDist <- dplyr::mutate(.data = rateDist,
                                  proportion = proportion * 100)
        
        estPar <- mdl[[m]][['estimated parameters']]
        
        tibble::tibble(
          'AIC-c' = aicc,
          'Log-likelihood' = loglike,
          'Rate Distribution' = list(rateDist),
          'Estimated parameters' = estPar
        )
      } else if (m == 'MG94xREV with separate rates for branch sets') {
        
        # Get Alpha-sorted codons (excluding STOP codons)
        codons <- sort(names(Biostrings::getGeneticCode()))
        codons <- codons[ ! codons %in% c('TAG', 'TAA', 'TGA') ]
        
        eqFreq <- magrittr::set_names(x = unlist(mdl[[m]][['Equilibrium frequencies']]),
                                      value = codons)
        
        aicc <- mdl[[m]][['AIC-c']]
        loglike <- mdl[[m]][['Log Likelihood']]
        rateDist <- unlist(mdl[[m]][['Rate Distributions']])
        estPar <- mdl[[m]][['estimated parameters']]
        
        tibble::tibble(
          'AIC-c' = aicc,
          'Equilibrium Frequency' = list(eqFreq),
          'Log-likelihood' = loglike,
          'Rate Distribution' = list(rateDist),
          'Estimated parameters' = estPar
        )
      } else if(m == 'Nucleotide GTR') {
        aicc <- mdl[[m]][['AIC-c']]
        eqFreq <-magrittr::set_names(x = unlist(mdl[[m]][['Equilibrium frequencies']]),
                                     value = c('A', 'C', 'G', 'T'))
        loglike <- mdl[[m]][['Log Likelihood']]
        rateDist <- unlist(mdl[[m]][['Rate Distributions']])
        estPar <- mdl[[m]][['estimated parameters']]
        
        tibble::tibble(
          'AIC-c' = aicc,
          'Equilibrium Frequency' = list(eqFreq),
          'Log-likelihood' = loglike,
          'Rate Distribution' = list(rateDist),
          'Estimated parameters' = estPar
        )
      } else {
        aicc <- mdl[[m]][['AIC-c']]
        loglike <- mdl[[m]][['Log Likelihood']]
        
        rateDist <- unlist(mdl[[m]][['Rate Distributions']])
        rateDist <- tibble::tibble(Names = names(rateDist), 
                                   Values = rateDist)
        rateDist <- tidyr::separate(data = rateDist,
                                    col = Names, 
                                    into = c('Branch set',
                                             'Omega rate classes',
                                             'propRate'), 
                                    sep = '\\.')
        rateDist <- dplyr::filter(.data = rateDist,
                                  `Branch set` == 'Test')
        rateDist <- tidyr::pivot_wider(data = rateDist,
                                       names_from = propRate,
                                       values_from = Values)
        rateDist <- dplyr::mutate(.data = rateDist,
                                  proportion = proportion * 100)
        
        estPar <- mdl[[m]][['estimated parameters']]
        
        tibble::tibble(
          'AIC-c' = aicc,
          'Log-likelihood' = loglike,
          'Rate Distribution' = list(rateDist),
          'Estimated parameters' = estPar
        )
      }
    })
    
    # Bind results of each model per sequence
    names(df_models) <- nms
    df_models <- dplyr::bind_rows(df_models, .id = 'Model')
    df_models <- dplyr::arrange(.data = df_models,
                                dplyr::desc(`Log-likelihood`))
  })
  
  # Return dataframe of all sequences run
  return(dplyr::bind_rows(out, .id = 'ID'))
}

# Get branch-attributes data
getBranchAttributesBusted <- function(l) {
  out <- purrr::map(l, function(json) {
    
    # Drop 'attributes' key from partitions (not needed)
    partitions <- json[['branch attributes']]
    partitions <- partitions[-length(partitions)]
    
    # Iterate over partitions
    df_partition <- purrr::map(partitions, function(pt) {
      
      # Iterate over species within partition
      df_species <- purrr::map(pt, function(sp) {
        
        nm <- names(sp)
        
        if ('constrained' %in% nm) {
          tibble::tibble(
            'MG94xREV with separate rates for branch sets' = sp[['MG94xREV with separate rates for branch sets']],
            'Nucleotide GTR' = sp[['Nucleotide GTR']],
            'Constrained' = sp[['constrained']],
            'Unconstrained' = sp[['unconstrained']]
          )
        } else {
          tibble::tibble(
            'MG94xREV with separate rates for branch sets' = sp[['MG94xREV with separate rates for branch sets']],
            'Nucleotide GTR' = sp[['Nucleotide GTR']],
            'Unconstrained' = sp[['unconstrained']]
          )
        }
      })
      
      dplyr::bind_rows(df_species, .id = 'Species')
    })
    
    dplyr::bind_rows(df_partition, .id = 'Partition')
  })
  
  return(dplyr::bind_rows(out, .id = 'ID'))
}

# Get Site log-likelihood information
getSiteLogLikelihoodBusted <- function(l) {
  out <- purrr::map(l, function(json) {
    models <- names(json[['Site Log Likelihood']])
    
    # Iterate over models
    df_logLike <- purrr::map(models, function(m) {
      tibble::tibble(
        'Model' = m,
        'Position' = seq(length(json[['Site Log Likelihood']][[m]][[1]])),
        'logLikelihood' = unlist(json[['Site Log Likelihood']][[m]][[1]])
      )
    })
    df_logLike <- dplyr::bind_rows(df_logLike)
    df_logLike <- tidyr::pivot_wider(data = df_logLike,
                                     names_from = Model,
                                     values_from = logLikelihood)
  })
  
  # Bind dataframes for each sequence
  dplyr::bind_rows(out, .id = 'ID')
}

# Get evidence ratio data
getEvidenceRatiosBusted <- function(l) {
  out <- purrr::map(l, function(json) {
    models <- names(json[['Evidence Ratios']])
    
    # Iterate over models - Will be empty if 'constrained' model NOT fit
    if(!purrr::is_empty(models)) {
      df_erat <- purrr::map(models, function(m) {
        tibble::tibble(
          'Model' = m,
          'Position' = seq(length(json[['Evidence Ratios']][[m]][[1]])),
          'Evidence ratio' = unlist(json[['Evidence Ratios']][[m]][[1]])
        )
      })
    } else {
      df_erat <- tibble::tibble(
        'Model' = c('constrained', 'optimized null'),
        'Position' = NA_integer_,
        'Evidence ratio' = NA_integer_
      )
    }
    
    df_erat <- dplyr::bind_rows(df_erat)
    df_erat <- tidyr::pivot_wider(data = df_erat,
                                  names_from = Model,
                                  values_from = `Evidence ratio`)
    
  })
  
  # Bind dataframes for each sequence
  dplyr::bind_rows(out, .id = 'ID')
}

# ---------------------------------------------------------------------------- #
# MEME Specific

getMleMeme <- function(l, workers = 4) {
  colNames <- c('alpha', 'beta-', 'p-',
                'beta+', 'p+', 'LRT',
                'p-value', '# branch selection',
                'total branch length', 'MEME Logl',
                'FEL Logl')
  
  future::plan('multisession', workers = workers)
  # TIME: 17.502
  lst_json <- furrr::future_map(l, function(json) {
    mle <- json[['MLE']][['content']]
    lst <- lapply(mle, function(partition) {
      d <- do.call(rbind, partition)
      colnames(d) <- colNames
      d <- tidyr::unnest(tibble::as_tibble(d), cols = dplyr::all_of(colNames))
    })
    dplyr::bind_rows(lst, .id = 'Partition')
  })
  return(dplyr::bind_rows(lst_json, .id = 'ID'))
}

# TIME: 2.765
getFitsMeme <- function(l) {
  out <- purrr::map(l, function(json) {
    
    # Store model names and models into objects
    nms <- names(json[['fits']])
    mdl <- json[['fits']]
    
    # Iterate over each model and extract relevant information
    df_models <- purrr::map(nms, function(m) {
      if (m == 'Global MG94xREV') {
        
        # Get Alpha-sorted codons (excluding STOP codons)
        codons <- sort(names(Biostrings::getGeneticCode()))
        codons <- codons[ ! codons %in% c('TAG', 'TAA', 'TGA') ]
        
        eqFreq <- magrittr::set_names(x = unlist(mdl[[m]][['Equilibrium frequencies']]),
                                      value = codons)
        
        aicc <- mdl[[m]][['AIC-c']]
        loglike <- mdl[[m]][['Log Likelihood']]
        rateDist <- unlist(mdl[[m]][['Rate Distributions']])
        estPar <- mdl[[m]][['estimated parameters']]
        
        tibble::tibble(
          'AIC-c' = aicc,
          'Equilibrium Frequency' = list(eqFreq),
          'Log-likelihood' = loglike,
          'Rate Distribution' = list(rateDist),
          'Estimated parameters' = estPar
        )
      } else if(m == 'Nucleotide GTR') {
        aicc <- mdl[[m]][['AIC-c']]
        eqFreq <-magrittr::set_names(x = unlist(mdl[[m]][['Equilibrium frequencies']]),
                                     value = c('A', 'C', 'G', 'T'))
        loglike <- mdl[[m]][['Log Likelihood']]
        rateDist <- unlist(mdl[[m]][['Rate Distributions']])
        estPar <- mdl[[m]][['estimated parameters']]
        
        tibble::tibble(
          'AIC-c' = aicc,
          'Equilibrium Frequency' = list(eqFreq),
          'Log-likelihood' = loglike,
          'Rate Distribution' = list(rateDist),
          'Estimated parameters' = estPar
        )
      }
    })
    
    # Bind results of each model per sequence
    names(df_models) <- nms
    df_models <- dplyr::bind_rows(df_models, .id = 'Model')
    df_models <- dplyr::arrange(.data = df_models,
                                dplyr::desc(`Log-likelihood`))
  })
  
  # Return dataframe of all sequences run
  return(dplyr::bind_rows(out, .id = 'ID'))
}

getBranchAttributesMeme <- function(l) {
  out <- purrr::map(l, function(json) {
    
    # Drop 'attributes' key from partitions (not needed)
    partitions <- json[['branch attributes']]
    partitions <- partitions[-length(partitions)]
    
    df_partition <- purrr::map(partitions, function(pt) {
      
      df_sp <- purrr::map(pt, function(sp) {
        
        l <- length(sp)
        ebf <- sp[1:(l - 2)]
        mdl <- sp[(l - 1):l]
        
        # Empirical Bayes Factor
        ebf <- tibble::tibble('Position' = names(ebf),
                              'Empirical Bayes Factor' = unlist(ebf))
        mdl <- tibble::tibble('Model' = names(mdl),
                              'Branch length' = unlist(mdl))
        
        tibble::tribble(
          ~'EBF',~'Model',
          ebf,mdl
        )
      })
      
      # Tibble of all species in analysis
      dplyr::bind_rows(df_sp, .id = 'Species')
    })
    
    # Tibble of all partitions in an analysis
    dplyr::bind_rows(df_partition, .id = 'Partition')
  })
  
  # Tibble across all sequences analysed
  return(dplyr::bind_rows(out, .id = 'ID'))
}

# ---------------------------------------------------------------------------- #
# Relax

getTestResultsRelax <- function(l) {
  out <- purrr::map(l, function(json) {
    tibble::tibble(
      'LRT' = json[['test results']][['LRT']],
      'p-value' = json[['test results']][['p-value']],
      'relaxation or intensification parameter' =json[['test results']][['relaxation or intensification parameter']]
    )
  })
  
  return(dplyr::bind_rows(out, .id = 'ID'))
}

getFitsRelax <- function(l) {
  out <- purrr::map(l, function(json) {
    df_models <- purrr::imap(json[['fits']], ~{
      
      # Equilibrium frequencies data
      codons <- sort(names(Biostrings::getGeneticCode()))
      codons <- codons[ ! codons %in% c('TAG', 'TAA', 'TGA') ]
      
      nuc <- c('A', 'C', 'G', 'T')
      
      # Common to all models
      aicc <- .x[['AIC-c']]
      logl <- .x[['Log Likelihood']]
      estPar <- .x[['estimated parameters']]
      rateDist <- .x[['Rate Distributions']]
      
      if (.y %in% c('General descriptive', 'RELAX alternative', 
                    'RELAX null', 'RELAX partitioned descriptive')) {
        tibble::tibble(
          'AIC-c' = aicc,
          'Log-likelihood' = logl,
          'Rate Distribution' = list(rateDist),
          'Estimated parameters' = estPar
        )
      } else if (.y == 'MG94xREV with separate rates for branch sets') {
        eqFreq <- magrittr::set_names(x = unlist(.x[['Equilibrium frequencies']]),
                                      value = codons)
        tibble::tibble(
          'AIC-c' = aicc,
          'Equilibrium Frequencies' = list(eqFreq),
          'Log-likelihood' = logl,
          'Rate Distribution' = list(rateDist),
          'Estimated parameters' = estPar
        )
        
      } else if(.y == 'Nucleotide GTR') {
        
        eqFreq <- magrittr::set_names(x = unlist(.x[['Equilibrium frequencies']]),
                                      value = nuc)
        tibble::tibble(
          'AIC-c' = aicc,
          'Equilibrium Frequencies' = list(eqFreq),
          'Log-likelihood' = logl,
          'Rate Distribution' = list(rateDist),
          'Estimated parameters' = estPar
        )
        
      }
    })
    
    dplyr::bind_rows(df_models, .id = 'Model')
    
  })
  
  return(dplyr::bind_rows(out, .id = 'ID'))
}

getBranchAttributesRelax <- function(l) {
  out <- purrr::map(l, function(json) {
    partitions <- json[['branch attributes']]
    partitions <- partitions[-length(partitions)]
    
    df_partition <- purrr::map(partitions, function(pt) {
      
      df_sp <- purrr::imap(pt, ~{
        tib <- tibble::tibble(
          columns = names(.x),
          values = unlist(.x)
        )
        
        tib <- dplyr::mutate(
          tib,
          Species = .y
        )
        
        tib <- tidyr::pivot_wider(
          tib,
          names_from = columns,
          values_from = values
        )
        tib <- dplyr::mutate(
          tib,
          dplyr::across(
            .cols = 2:8,
            .fns = as.numeric
          )
        )
      })
      
      df_sp <- dplyr::bind_rows(df_sp)
      df_sp <- dplyr::select(df_sp, -`original name`)
    })
    
    dplyr::bind_rows(df_partition, .id = 'Parition')
    
  })
  
  return(dplyr::bind_rows(out, .id = 'ID'))
}
