################################################################################
# This is a temporary script that I've pushed to GitHub in order to have a copy
# in version history. It compares the probablities computed using the
# approximate joint before and after setting recrudescence probabilities to zero
# for non-adjacent episodes — a modification following discussion with
# YSF on 1st April 2026.
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

# These are identical: 
identical(unlist(Uniform_pairwise_old), unlist(Uniform_pairwise))

# Because the posterior is either zero or greater than 0.5 and
# compute_posterior_approxjoint contains the if clause for recrudescence probs
# greater than 0.5
C_probs <- unlist(sapply(Uniform_pairwise, function(x) x[,"C"]))
range(C_probs[C_probs > 0.5])
range(C_probs[C_probs <= 0.5])

# These are not: 
identical(unlist(TimeToEvent_pairwise_old), unlist(TimeToEvent_pairwise))
max(abs(unlist(TimeToEvent_pairwise_old) - unlist(TimeToEvent_pairwise)))
C_probs <- unlist(sapply(TimeToEvent_pairwise, function(x) x[,"C"]))
any(C_probs > 0.5) # No C probs larger than 0.5
range(C_probs[C_probs <= 0.5]) # Contains some non zero C probs

