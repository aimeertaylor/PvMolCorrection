################################################################################
# Allele frequency estimates
#
# Our approach undercounts repeat alleles in episodes with MOIs > 1 because it
# ignores NAs and NAs are used as fillers at markers with observed cardinality <
# MOI (see Understanding_NAs.R)
#
# In the Taylor & Watson et al. 2019 we intended to estimate population-level
# frequencies using enrolment episodes only but unintentionally used all
# episodes (because Ind_Primary was not passed to apply). Given most recurrences
# are likely to be either reinfections or relapses, both of which are draws from
# the mosquito population (albeit a delayed draw in the case of a relapse),
# in the absence of systematic within-patient selection (as might occur when
# break-through infections encounter lingering drug pressure), estimates based
# on all episodes should be unbiased and more precise than those based on
# enrolment episodes only. As such, it is not necessarily a problem that we
# didn't fulfil our original intention to exclude data from recurrent
# infections.
################################################################################
rm(list = ls())
library(Pv3Rs) # To get fs_VHX_BPD

# Old allele frequencies: FS_combined
load('~/Documents/RecurrentVivax/RData/Data_for_relatedness.RData')

# Compare old and fs_VHX_BPD: doesn't appear to be any systematic bias consistent with
# most samples being reinfections/relapses not subject to any systematic
# within-patient selection
for(MS in names(fs_VHX_BPD)){
  plot(fs_VHX_BPD[[MS]], Fs_Combined[[MS]])
  abline(a = 0, b = 1, lty = 'dashed')
}

