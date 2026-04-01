################################################################################
# Check number of markers typed
################################################################################
library(Pv3Rs)
rm(list = ls())
source("match_counting.R")

# Tag participants with more than one episode typed
rep_epi_log <- sapply(ys_VHX_BPD, length) > 1 

# Partition markers into those that were genotyped first and then extra
markers <- names(fs_VHX_BPD)
markers_3 <- c("PV.3.27", "PV.3.502", "PV.ms8") 
markers_4plus <- setdiff(markers, markers_3)

# PIDs with only three markers typed across two or more episodes
pids_3_log <- sapply(ys_VHX_BPD, function(y) {
  all(sapply(y, function(y_epi) !names(y_epi) %in% markers_4plus)) # No extra markers typed
 })
pids_3 <- names(which(pids_3_log & rep_epi_log)) 
matches_3 <- sapply(ys_VHX_BPD[pids_3], match_counting) # Compute allele matches_3
matches_3_p <- unlist(sapply(matches_3, function(x) x["p",])) # Extract proportion for all recurrences
table(matches_3_p) # Majority have <= 1/3 matches_3:
sum(matches_3_p <= 1/3)/length(matches_3_p) # 88 %
excep <- names(which(matches_3_p  > 1/3)) # exceptions to <= 1/3
pids_3_excep <- unique(do.call(rbind, strsplit(excep, split = ".", fixed = T))[,1])
plot_data(ys_VHX_BPD[pids_3_excep], fs = fs_VHX_BPD)  
matches_3[pids_3_excep] # fewer than two match between 1st and 2nd episode

# Is the rule: if fewer than two match between the first and second episode,
# type only three markers; otherwise, attempt all markers? Not quite, some
# participants have data on more than three makers despite fewer than two
# matches between the first and second episode... 
pids_4plus <- setdiff(names(ys_VHX_BPD)[rep_epi_log], pids_3)
plot_data(ys_VHX_BPD[pids_4plus], fs = fs_VHX_BPD)  

# Compute match proportion among key markers for pids_4plus
matches_4plus_p <- unlist(sapply(ys_VHX_BPD[pids_4plus], function(y) {
  z <- sapply(y, function(y_epi) y_epi[names(y_epi)[names(y_epi) %in% markers_3]], simplify = F)
  if(all(sapply(z, length)[1:2] > 0)) { # If data on at least one key marker across episodes 1 and 2
    p <- match_counting(z)["p",1] # This is only the first recurrence
  } else {
    p <- NA
  }})) 

# Participants where key markers failed for first episode: 
plot_data(ys_VHX_BPD[names(which(is.na(matches_4plus_p)))], fs = fs_VHX_BPD)

# Participants that have data on more than three markers despite fewer than two
# matches between the first and second episode...
pids_4plus_excep <- names(which(matches_4plus_p[!is.na(matches_4plus_p)] <= 1/3))
plot_data(ys_VHX_BPD[pids_4plus_excep], fs = fs_VHX_BPD)

