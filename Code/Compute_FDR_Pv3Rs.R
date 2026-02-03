rm(list = ls())
library(Pv3Rs) # For plot_data
library(tictoc)
Run_data <- FALSE # Run data formating versus loading pre-formatted data
Run_probs <- FALSE # Run posterior probability computation (2 hrs) versus loading pre-computed results
Run_irred_prob <- FALSE # Run irreducible probability computation 
tot_MOI_lim = 8 # Upper limit on computation
n_mark_lim = 2 # Lower limit on posterior computation
#===============================================================================
# Format data for Pv3Rs
#===============================================================================
if (!Run_data) {
  load("../RData/ys_FRD.RData")
} else { # Format APC_MSdata into Pv3Rs input:
  load('~/Documents/RecurrentVivax/RData/LargeFiles/APC_MSdata.bigRData')
  MSs_all <- colnames(APC_MSdata)[grepl("PV.", colnames(APC_MSdata))] # Extract marker names
  y.dfs <- plyr::dlply(APC_MSdata, 'ID') # Group genetic data by patient ID
  names(y.dfs) <- gsub("%EP2", "", gsub("APC_EP1%", "", names(y.dfs))) # Rename
  n_test <- 100
  tic()
  ys_FDR <- lapply(names(y.dfs), function(pid) {
    
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
    
    # Covert episode IDs into episode indices
    names(y) <- c(1,2)
    
    # Return patient data
    return(y)
  })
  toc()  
  save(ys_FDR, file = "../RData/ys_FRD.RData")
}

#===============================================================================
# Run Pv3Rs
#===============================================================================
if (!Run_probs){
  load("../RData/probs_FDR.RData")
} else {
  n_test <- length(ys_FDR) # 250303 comparisons in ys_FDR
  tic()
  probs <- sapply(1:n_test, function(i){
    tot_MOI <- sum(Pv3Rs::determine_MOIs(ys_FDR[[i]]))
    n_mark <- length(intersect(names(ys_FDR[[i]][[1]]), names(ys_FDR[[i]][[2]])))
    if(tot_MOI <= tot_MOI_lim & n_mark > n_mark_lim) {
      return(suppressMessages(Pv3Rs::compute_posterior(y = ys_FDR[[i]], fs = fs_VHX_BPD)$marg))
    } else {
      return(rep(NA,3)) 
    }
  })
  toc()
  save(probs, file = "../RData/probs_FDR.RData")
}


#===============================================================================
# Run Pv3Rs
#===============================================================================
if (!Run_irred_prob){
  load("../RData/irred_probs_FDR.RData")
} else {
  load("../../Pv3Rs/vignettes/articles/MaximaStudy/maxima.rda") 
  n_test <- length(ys_FDR) # 250303 comparisons in ys_FDR
  tic()
  irred_probs <- sapply(1:n_test, function(i){
    MOIs <- Pv3Rs::determine_MOIs(ys_FDR[[i]])
    if(sum(MOIs) <= tot_MOI_lim) {
      return(1-maxima["I_with", paste(MOIs, collapse = "")])
    } else {
      return(NA)
    }
  })
  toc()
  save(irred_probs, file = "../RData/irred_probs_FDR.RData")
}


#===============================================================================
# Compute false failure discovery rate: 
#===============================================================================
uncomputable <- is.na(probs[3,])
Pv3Rs::plot_data(ys_FDR[uncomputable], fs = fs_VHX_BPD) # Only one marker with data
n_computable <- sum(!uncomputable) # Denominator for false discovery rates 

# Irreducible Failure probability rate
100 * round(sum(irred_probs[!uncomputable])/n_computable, 3)

# Failure probability rate:
100 * round(sum(1-probs[3,], na.rm = T)/n_computable, 3)

# Failure classification rate; see lines 56 and 1656 of Pooled_Analysis.Rmd
eps <- 0.7 # O.7 on line 56 of Pooled_Analysis.Rmd
100 * round(sum((1-probs[3,]) > eps, na.rm = T)/n_computable, 3)

