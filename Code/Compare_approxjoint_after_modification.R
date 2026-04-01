################################################################################
# This is a temporary script that I've pushed to GitHub in order to have a copy
# in version history. It compares the probabilities computed using the
# approximate joint before and after 
# 1) setting recrudescence probabilities to zero
# 2) implementing Elijah's change suggested 1st April 2026
################################################################################
rm(list = ls())

load("../RData/results_Pv3Rs_old.RData")
TimeToEvent_pairwise_old <- TimeToEvent_pairwise
Uniform_pairwise_old <- Uniform_pairwise
load("../RData/results_Pv3Rs.RData")

# Comparing results after setting recrudescent priors to zero
par(mfrow = c(2,1))
plot(unlist(TimeToEvent_pairwise_old), unlist(TimeToEvent_pairwise))
plot(unlist(Uniform_pairwise_old), unlist(Uniform_pairwise))

# These are no-longer identical: 
identical(unlist(Uniform_pairwise_old), unlist(Uniform_pairwise))
max(abs(unlist(Uniform_pairwise_old) - unlist(Uniform_pairwise)))
which(abs(unlist(Uniform_pairwise_old) - unlist(Uniform_pairwise)) > 0.01)
Uniform_pairwise_old[["VHX_220"]]
Uniform_pairwise[["VHX_220"]]
plot_data(y)

# These are not: 
identical(unlist(TimeToEvent_pairwise_old), unlist(TimeToEvent_pairwise))
max(abs(unlist(TimeToEvent_pairwise_old) - unlist(TimeToEvent_pairwise)))
which(abs(unlist(TimeToEvent_pairwise_old) - unlist(TimeToEvent_pairwise)) > 0.01)
plot_data(ys_VHX_BPD[c("VHX_350", "VHX_646")])
TimeToEvent_pairwise_old["VHX_646"]
TimeToEvent_pairwise["VHX_646"]

