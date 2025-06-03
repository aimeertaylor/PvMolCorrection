################################################################################
# Role of NAs in the microsatellite data: NAs encode missingness and are used as
# fillers for markers with allele counts below the MOI
################################################################################
rm(list = ls())

# Genetic data (MS_pooled)
load('../jwatowatson-RecurrentVivax-4870715/RData/GeneticModel/MS_data_PooledAnalysis.RData') 

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

# Min allele count zero: confirms NAs encode missingness
min(within_episode_allele_counts)

# Max allele count one (no repeats): confirms NAs are used as filler
max(within_episode_allele_counts)

# Example of NA use as filler: 
MS_pooled_byepisode$VHX_225_3
