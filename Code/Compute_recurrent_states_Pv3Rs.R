################################################################################
# In this script estimates are generated under the Pv3Rs model using:
#
# 1) a uniform prior
# 2) normalised time-to-event posterior mean estimates
#
# To compute median-based PQ failure rates with CIs need to generate estimates
# for all 100 posterior time-to-event draws stored in
# load('../RData/TimingModel/MOD2_Posterior_samples.RData') and 100 draws from
# the posterior allele frequency distribution, which we cannot replicate exactly
# because we didn't use set.seed in Pooled_Analysis.Rmd. This would cost approx.
# 100 hours of computation - an unjustified expense given the focus here is on 
# the genetic analysis.
################################################################################

rm(list = ls())
library(Pv3Rs)
source("compute_posterior_approxjoint.R")

#===============================================================================
# Get prior estimates based on time-to-event and normalise
#===============================================================================
# Load time-to-event estimates generated in original study 
load("../jwatowatson-RecurrentVivax-4870715/RData/TimingModel/MOD2_theta_estimates.RData")

# Extract posterior means (Pooled_Analysis.Rmd lines 672 and 1807)
inds = grepl('mean_theta', colnames(Mod2_ThetaEstimates))
prior_unnorm = data.frame(Mod2_ThetaEstimates[,inds], stringsAsFactors = F)
rownames(prior_unnorm) <- Mod2_ThetaEstimates$Episode_Identifier

# Rename columns
colnames(prior_unnorm) = gsub(pattern = 'Recrudescence_mean_theta', replacement = 'C', x = colnames(prior_unnorm))
colnames(prior_unnorm) = gsub(pattern = 'Relapse_mean_theta', replacement = 'L', x = colnames(prior_unnorm))
colnames(prior_unnorm) = gsub(pattern = 'ReInfection_mean_theta', replacement = 'I',x = colnames(prior_unnorm))

# Normalise s.t. all prior estimates sum to one (this step was skipped before)
prior <- prior_unnorm / rowSums(prior_unnorm)

# Save prior estimates for DevFiles/Understanding_prior.R
save(prior, prior_unnorm, file = "../RData/prior_estimates.RData")

# ==============================================================================
# Compute new results
# ==============================================================================
TimeToEvent_joint <- Uniform_joint <- TimeToEvent_pairwise <- Uniform_pairwise <- list()
paired_data_pids <- names(which(sapply(ys_VHX_BPD, length) > 1))
tictoc::tic()
for(pid in paired_data_pids) {
  
  writeLines(pid)
  
  # Extract data and prior for recurrent episodes
  y <- ys_VHX_BPD[[pid]]
  sum_MOIs <- sum(determine_MOIs(y))
  prior_per_patient <- prior[paste0(pid, "_", names(y)[-1]),]
  row.names(prior_per_patient) <- names(y)[-1] # Re-name episodes in prior 
  
  # Store result in a list
  if (sum_MOIs > 8) { # Don't compute joint if more than 8 genotypes
    TimeToEvent_joint[[pid]] <- NA
    Uniform_joint[[pid]] <- NA
    TimeToEvent_pairwise[[pid]] <- compute_posterior_approxjoint(y, fs_VHX_BPD, prior = prior_per_patient)
    Uniform_pairwise[[pid]] <- compute_posterior_approxjoint(y, fs_VHX_BPD)    
  } else { 
    TimeToEvent_joint[[pid]] <- compute_posterior(y, fs_VHX_BPD, prior = prior_per_patient)
    Uniform_joint[[pid]] <- compute_posterior(y, fs_VHX_BPD)
    TimeToEvent_pairwise[[pid]] <- compute_posterior_approxjoint(y, fs_VHX_BPD, prior = prior_per_patient)
    Uniform_pairwise[[pid]] <- compute_posterior_approxjoint(y, fs_VHX_BPD)
  }
}
tictoc::toc()

# Save
save(TimeToEvent_joint, Uniform_joint, TimeToEvent_pairwise, Uniform_pairwise,
     file = "../RData/results_Pv3Rs.RData")
