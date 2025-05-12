################################################################################
# Confirming the role of NAs in the microsatellite data: Intra-episode per
# marker allele counts of zero confirm NAs encode missingness. Intra-episode per
# marker maximum allele counts of one confirm NAs are also used as fillers for
# markers whose observed cardinality is less than the MOI when MOI > 1
################################################################################
rm(list = ls())

# Genetic data (MS_pooled)
load('~/Documents/RecurrentVivax/RData/GeneticModel/MS_data_PooledAnalysis.RData') 

# Marker names
MSs_all <- tail(names(MS_pooled), 9) 

# Split data into episodes 
MS_pooled_byepisode <- plyr::dlply(MS_pooled, "MS_pooled$Episode_Identifier")

# For each episode, for each marker, tabulate allele observations 
within_episode_allele_counts <- sapply(MS_pooled_byepisode, function(x) {
  list_allele_counts <- apply(x[MSs_all], 2, table)
  per_marker_allele_counts <- sapply(list_allele_counts, function(count) {
    if (length(count) == 0) {
      0
    } else {
      max(count)
    }
  })
})

# Max allele count is one (confirms NAs are used as a filler)
range(within_episode_allele_counts)

# Example 
MS_pooled_example <- MS_pooled_byepisode$VHX_225_3
