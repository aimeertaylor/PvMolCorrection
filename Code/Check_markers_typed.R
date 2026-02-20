################################################################################
# Check number of markers typed
################################################################################
rm(list = ls())
source("match_counting.R")

# PIDs with only three markers successfully typed across two or more episodes
pids_3_log <- unlist(sapply(ys_VHX_BPD, function(y) all(sapply(y, length) == 3)))
rep_epi_log <- sapply(ys_VHX_BPD, length) > 1
pids <- names(which(pids_3_log & rep_epi_log))
matches <- sapply(ys_VHX_BPD[pids], match_counting) # Compute allele matches
matches_p <- unlist(sapply(matches, function(x) x["p",])) # Extract proportion
table(matches_p) # Majority have <= 1/3 matches:
sum(matches_p <= 1/3)/length(matches_p) # 88

# Exceptions: 
excep <- names(which(matches_p  > 1/3))
pids_excep <- unique(do.call(rbind, strsplit(excep, split = ".", fixed = T))[,1])
all(pids_excep %in% names(ys_VHX_BPD))
plot_data(ys_VHX_BPD[pids_excep], fs = fs_VHX_BPD)  
