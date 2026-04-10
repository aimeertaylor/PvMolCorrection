rm(list = ls())
load("../RData/results_Pv3Rs.RData")
par(mfrow = c(1,2))

# Extract pids for people with jointly inferred probabilities: 
# same participants, regardless of prior
pids_joint <- names(TimeToEvent_joint)[!is.na(TimeToEvent_joint)]
identical(pids_joint, names(Uniform_joint)[!is.na(Uniform_joint)])

# Find participants with more than one recurrence
num_epi <- sapply(ys_VHX_BPD[pids_joint], length) 
pids_joint_seq <- names(num_epi)[num_epi > 2]
length(pids_joint_seq) # 57 people

# ==============================================================================
# Uniform:
# ==============================================================================
# Extract probability of failure for each: 
failure_prob <- sapply(Uniform_joint[pids_joint_seq], function(x) {
  i <- which(grepl("C", names(x$joint)) + grepl("L", names(x$joint)) == 0)
  c(indep_not_assum = unname(1-x$joint[i]), 
    indep_assum = 1-prod(x$marg[,"I"]))
})

plot(x = failure_prob["indep_not_assum",], 
     y = failure_prob["indep_assum",], bty = "n", pch = 20, 
     xlab = "Independence not assumed", ylab = "Independence assumed", 
     main = "Uniform prior")
abline(a = 0, b = 1, lty = "dashed")

max(abs(failure_prob["indep_not_assum",] - failure_prob["indep_assum",]))
sum(abs(failure_prob["indep_not_assum",] - failure_prob["indep_assum",]))

# ==============================================================================
# Time-to-event: 
# ==============================================================================
# Extract probability of failure for each: 
failure_prob <- sapply(TimeToEvent_joint[pids_joint_seq], function(x) {
  i <- which(grepl("C", names(x$joint)) + grepl("L", names(x$joint)) == 0)
  c(indep_not_assum = unname(1-x$joint[i]), 
    indep_assum = 1-prod(x$marg[,"I"]))
})

plot(x = failure_prob["indep_not_assum",], 
     y = failure_prob["indep_assum",], bty = "n", pch = 20, 
     xlab = "Independence not assumed", ylab = "Independence assumed", 
     main = "Time-to-event prior")
abline(a = 0, b = 1, lty = "dashed")

max(abs(failure_prob["indep_not_assum",] - failure_prob["indep_assum",]))
sum(abs(failure_prob["indep_not_assum",] - failure_prob["indep_assum",]))