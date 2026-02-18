# Plot per-marker allele count distribution ordered by marker cardinality
rm(list = ls())
library(Pv3Rs)

# Number of alleles per marker per participant
markers <- names(fs_VHX_BPD)
# Make store for allele counts per marker
allele_count_per_marker <- rep(NA, length(markers)); names(allele_count_per_marker) <- markers                               
# Extract allele counts per marker
allele_counts <- sapply(ys_VHX_BPD, function(y) {
  x <- sapply(y[[1]], length)
  allele_count_per_marker[names(x)] <- x
  allele_count_per_marker
})

# Marker cardinality:
eff_card <- sapply(fs_VHX_BPD, function(m) 1/sum(m^2))
# Make store for allele count categorical distribution:
allele_count_dist <- rep(NA,4); names(allele_count_dist) <- 1:4
# Extract allele count categorical distribution
allele_dist <- apply(allele_counts, 1, function(x){
  z <- table(x)
  allele_count_dist[names(z)] <- z
  allele_count_dist
})
# Convert to allele count frequency
allele_freq <- apply(allele_dist, 2, function(x) x/sum(x, na.rm = T))

# Plot distribution: no obvious outliers
barplot(allele_freq[,names(sort(eff_card))], beside = T)
