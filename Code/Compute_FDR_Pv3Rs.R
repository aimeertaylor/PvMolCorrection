rm(list = ls())
library(Pv3Rs) # For plot_data
library(tictoc)
Run_data <- F # Run data formatting versus loading pre-formatted data
Run_probs <- F # Run posterior probability computation (2 hours) versus loading pre-computed results
Run_min_prob <- F # Run minimum (irreducible) probability computation 
tot_MOI_lim = 8 # Upper limit on computation

#===============================================================================
# Format data for Pv3Rs
#===============================================================================
if (Run_data) {
  
  load('~/Documents/RecurrentVivax/RData/LargeFiles/APC_MSdata.bigRData')
  MSs_all <- colnames(APC_MSdata)[grepl("PV.", colnames(APC_MSdata))] # Extract marker names
  y.dfs <- plyr::dlply(APC_MSdata, 'ID') # Group genetic data by patient ID
  
  tic()
  ys_null <- lapply(names(y.dfs), function(pid) {
    
    # Extract data for a given patient IDs
    y.df <- y.dfs[[pid]]
    
    # Extract markers for which there is at least one non-NA
    ms <- MSs_all[apply(!is.na(y.df[MSs_all]), 2, any)]
    
    # Group data by episode
    y.by.episode <- plyr::dlply(y.df, 'Episode_Identifier')
    
    # Transform data frame format to format taken by 'compute_posterior'
    y <- lapply(y.by.episode, function(episode) { # For each episode
      setNames(lapply(ms, function(m) { # For each marker
        alleles <- episode[m][!is.na(episode[m])] # Extract non-NAs
        if (length(alleles) > 0) {
          return(alleles)
        } else {
          return(NA)
        }
      }), ms)
    })
    names(y) <- c(1,2) # Covert episode IDs into episode indices
    return(y) # Return patient data
  })
  toc()  
  
  # Name the null pairs and save
  names(ys_null) <- names(y.dfs)
  save(ys_null, file = "../RData/ys_null.RData")
  
} else {
  load("../RData/ys_null.RData")
} 

#===============================================================================
# Run Pv3Rs
#===============================================================================
if (Run_probs){
  tic()
  probs <- t(sapply(ys_null, function(y){
    tot_MOI <- sum(Pv3Rs::determine_MOIs(y))
    n_mark <- length(intersect(names(y[[1]]), names(y[[2]])))
    if(tot_MOI <= tot_MOI_lim) {
      return(suppressMessages(Pv3Rs::compute_posterior(y, fs = fs_VHX_BPD)$marg[1,]))
    } else {
      return(c("C" = NA, "L" = NA, "I" = NA)) 
    }
  }))
  toc()
  save(probs, file = "../RData/probs_null.RData")
} else {
  load("../RData/probs_null.RData")
} 


#===============================================================================
# Extract upper bounds on reinfection
#===============================================================================
if (Run_min_prob){
  load("../../Pv3Rs/vignettes/articles/MaximaStudy/maxima.rda") 
  tic()
  min_probs <- sapply(ys_null, function(y){
    MOIs <- Pv3Rs::determine_MOIs(y)
    if(sum(MOIs) <= tot_MOI_lim) {
      return(1-maxima["I_with", paste(MOIs, collapse = "")])
    } else {
      return(NA)
    }
  })
  toc()
  save(min_probs, file = "../RData/min_probs_null.RData")
} else {
  load("../RData/min_probs_null.RData")
} 


#===============================================================================
# Explore subsets of comparisons used in false discovery rate computations
#===============================================================================
# 250303 comparisons: 
n_comp <- length(ys_null) 

# For each pair compute number of markers with allele data on both episodes:
n_marks <- sapply(ys_null, function(y){
  n_mark <- sum(sapply(y[[1]], function(x) !all(is.na(x))) & 
                  sapply(y[[2]], function(x) !all(is.na(x))))  
})
one_plus <- names(which(n_marks > 0)) # not all comparisons have comparable data
n_one_plus <- length(one_plus) # 247873 comparisons with data one plus marker

# When no markers with common data, both Pv3Rs and the prototype return 
# probabilities close to but not equal to the prior; see link below 
# https://aimeertaylor.github.io/Pv3Rs/articles/posterior-probabilities.html#incomparable
load('~/Documents/RecurrentVivax/RData/LargeFiles/Inflated_Results.bigRData')
Inflated_Results[which(n_marks == 0)[1],]
suppressMessages(Pv3Rs::compute_posterior(ys_null[[which(n_marks == 0)[1]]], fs = fs_VHX_BPD)$marg)

# 249540 comparisons reported in Nat Comms (UpperComplexity = 10^6 line 1640 of Pooled_Analysis.Rmd):
# UpperComplexity = 10^6 approximately equal to MOItot = 5
ys_null_MOItot <- sapply(ys_null, function(x) sum(determine_MOIs(x)))
max(ys_null_MOItot) # At most eight
length(ys_null) - sum(ys_null_MOItot > 5) # 249155 

# Check numbers of comparisons
length(min_probs) == n_comp & sum(!is.na(probs[,1])) == n_comp
nrow(Inflated_Results) # smaller

#===============================================================================
# Compute false failure discovery rates with and without =comparisons with no
# data on common markers (makes little difference):
#===============================================================================
eps <- 0.7 # O.7 on line 56 of Pooled_Analysis.Rmd

# Irreducible (minimum) failure probability rate
100 * round(sum(min_probs)/n_comp, 3) # all comparisons
100 * round(sum(min_probs[one_plus])/n_one_plus, 3) # comparisons with data on one plus markers

# Failure probability rate:
100 * round(sum(1-probs[, "I"], na.rm = T)/n_comp, 3) # all comparisons
100 * round(sum(1-probs[one_plus, "I"], na.rm = T)/n_one_plus, 3) # only comparisons with data on one plus markers

# Failure classification rate using 2019 epsilon; see lines 56 and 1656 of Pooled_Analysis.Rmd
100 * round(sum((1-probs[, "I"]) > eps, na.rm = T)/n_comp, 3) # all comparisons
100 * round(sum((1-probs[one_plus, "I"]) > eps, na.rm = T)/n_one_plus, 3) # only comparisons with data on one plus markers

# Previous calculation
common_comps <- intersect(paste0(one_plus, "_2"), rownames(Inflated_Results))
100 * round(sum((1-Inflated_Results$I) > eps, na.rm = T)/sum(!is.na(1-Inflated_Results$I)), 3) # Reported 
100 * round(sum((1-Inflated_Results[common_comps, "I"]) > eps, na.rm = T)/sum(!is.na(Inflated_Results[common_comps, "I"])), 3) 


