################################################################################
# Allele frequency estimates
#
# This approach undercounts repeat alleles in episodes with MOIs > 1 because it
# ignores NAs and NAs are used as fillers at markers with observed cardinality <
# MOI (see Understanding_NAs.R)
#
# In the Taylor & Watson et al. 2019 we intended to estimate population-level
# frequencies using enrolment episodes only but unintentionally used all
# episodes (because Ind_Primary was not passed to apply). Given most recurrences
# are likely to be either reinfections or relapses, both of which are draws from
# the mosquito population (albeit a delayed draw in the case of a relapse),
# in the absence of systematic within-patient selection (as might occur when
# break through infections encounter lingering drug pressure), estimates based
# on all episodes should be unbiased and more precise than those based on
# enrolment episodes only. As such, it is not necessarily a problem that we
# didn't fulfil our original intention to exclude data from recurrent
# infections.
################################################################################
rm(list = ls())

# Genetic data (MS_pooled)
load('~/Documents/RecurrentVivax/RData/GeneticModel/MS_data_PooledAnalysis.RData') 

# Marker names
MSs_all <- tail(names(MS_pooled), 9) 

# Old allele frequencies (FS_combined)
load('../RData/Data_for_relatedness.RData') 

# Function to compute posterior counts
compute_posterior <- function(x, Indices){
  # number of possible alleles for this marker
  xmax <-  max(x, na.rm = T) 
  # prior concentration parameters
  param_vector <- array(D_weight_Prior, dim=xmax, dimnames=list(1:xmax))
  # observed data summarised as counts
  obs_counts <- table(x[Indices])
  # posterior concentration parameters
  param_vector[names(obs_counts)] <-  param_vector[names(obs_counts)] + obs_counts
  # end of function
  return(param_vector)
}

# Rename old 
Fs_Combined_old <- Fs_Combined[MSs_all] 

# Compute new from enrolment episodes only 
D_weight_Prior <- 1 # uniform Dirichlet prior for allele frequencies

# Recompute frequencies with enrolment episodes only 
Ind_enrol = which(MS_pooled$Episode==1) # Indices of primary episodes
Alpha_Posteriors <- apply(MS_pooled[,MSs_all], 2, compute_posterior, Ind_enrol)
Fs_Combined_enrol  <- sapply(Alpha_Posteriors, function(x){x/sum(x)})[MSs_all]

# Recompute frequencies with all episodes
Ind_all = rep(TRUE, nrow(MS_pooled)) # Indices of all episodes 
Alpha_Posteriors <- apply(MS_pooled[,MSs_all], 2, compute_posterior, Ind_all)
Fs_Combined_all  <- sapply(Alpha_Posteriors, function(x){x/sum(x)})[MSs_all]

# Compare old and new: doesn't appear to be any systematic bias consistent with
# most samples being reinfections/relapses not subject to any systematic
# within-patient selection
for(MS in MSs_all){
  plot(Fs_Combined_enrol[[MS]], Fs_Combined_old[[MS]])
  points(Fs_Combined_all[[MS]], Fs_Combined_old[[MS]], pch = 4)
  abline(a = 0, b = 1, lty = 'dashed')
}

