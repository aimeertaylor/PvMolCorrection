################################################################################
# Pairwise probabilities: lower probability of recrudescence; half-sib problem  
################################################################################

library(Pv3Rs)
load("../RData/results_Pv3Rs.RData")


# ==============================================================================
# Uniform: 
# Joint systematically has a slightly higher probability of recrudescence versus
# relapse. When the model is misspecified due to half siblings, a large
# difference between joint and pairwise estimates occurs if the model is not
# misspecified given one more episode pairs.
# ==============================================================================
# Extract marginal probabilities
Uniform_joint <- sapply(Uniform, function(x) x[[1]])

# Tag pids that were not analysed jointly
log_NA <- sapply(Uniform, function(x) any(is.na(x[[1]])))

# Remove pids that were not analysed jointly
Uniform_joint <- Uniform_joint[!log_NA]
Uniform_pwise <- Uniform_pairwise[!log_NA]

# Create vectors of estimates and features for plotting
unlst_uni_joint <- unlist(Uniform_joint)
unlst_uni_pwise <- unlist(Uniform_pwise)
n_rec <- sapply(Uniform_pwise, nrow) # Number of recurrences per pid
unlst_n_rec <- rep(rep(n_rec, n_rec), each = 3)
unlst_n_chr <- unlist(sapply(n_rec, function(x) rep(c("C","L","I"), each = x)))
unlst_names <- rep(rep(names(n_rec), n_rec), each = 3) 

# Plot discrepancy
plot(unlst_uni_joint, unlst_uni_pwise, ylab = "Pairwise", xlab = "Joint",
     col = unlst_n_rec, pch = unlst_n_chr, bty = "n")
legend("left", legend = unique(n_rec), fill = unique(n_rec), bty = "n", 
       title = "Recurrence \n count")
abline(a = 0.23, b = 1, lty = "dashed")
abline(a = -0.23, b = 1, lty = "dashed")
abline(a = 0, b = 1, lty = "dashed")

# Extract data and estimates for pid where joint and pairwise differ
big_diff_log <- abs(unlst_uni_joint - unlst_uni_pwise) > 0.23
pids_big_diff <- unique(unlst_names[big_diff_log])
prob_big_diff <- sapply(pids_big_diff, function(pid) {
  c(joint = unname(Uniform_joint[pid]), 
    pairwise = unname(Uniform_pwise[pid]))}, simplify = F)

# Inspect data and estimates where joint and pairwise differ
plot_data(ys = ys_VHX_BPD[pids_big_diff], fs = fs_VHX_BPD)
prob_big_diff


# ==============================================================================
# TimeToEvent: only one not-so-big (0.32) discrepency; lower relapse under
# joint, which seems more reasonable
# ==============================================================================
# Extract marginal probabilities
TimeToEvent_joint <- sapply(TimeToEvent, function(x) x[[1]])

# Tag pids that were not analysed jointly
log_NA <- sapply(TimeToEvent, function(x) any(is.na(x[[1]])))

# Remove pids that were not analysed jointly
TimeToEvent_joint <- TimeToEvent_joint[!log_NA]
TimeToEvent_pwise <- TimeToEvent_pairwise[!log_NA]

# Create vectors of estimates and features for plotting
unlst_uni_joint <- unlist(TimeToEvent_joint)
unlst_uni_pwise <- unlist(TimeToEvent_pwise)
n_rec <- sapply(TimeToEvent_pwise, nrow) # Number of recurrences per pid
unlst_n_rec <- rep(rep(n_rec, n_rec), each = 3)
unlst_n_chr <- unlist(sapply(n_rec, function(x) rep(c("C","L","I"), each = x)))
unlst_names <- rep(rep(names(n_rec), n_rec), each = 3) 

# Plot discrepancy
plot(unlst_uni_joint, unlst_uni_pwise, 
     col = unlst_n_rec, pch = unlst_n_chr, bty = "n", ylab = "Pairwise", xlab = "Joint")
legend("left", legend = unique(n_rec), fill = unique(n_rec), bty = "n", 
       title = "Recurrence \n count")
abline(a = 0.23, b = 1, lty = "dashed")
abline(a = -0.23, b = 1, lty = "dashed")
abline(a = 0, b = 1, lty = "dashed")

# Extract data and estimates for pid where joint and pairwise differ
big_diff_log <- abs(unlst_uni_joint - unlst_uni_pwise) > 0.23
pids_big_diff <- unique(unlst_names[big_diff_log])
prob_big_diff <- sapply(pids_big_diff, function(pid) {
  c(joint = unname(TimeToEvent_joint[pid]), 
    pairwise = unname(TimeToEvent_pwise[pid]))}, simplify = F)

# Inspect data and estimates where joint and pairwise differ
plot_data(ys = ys_VHX_BPD[pids_big_diff], fs = fs_VHX_BPD)
prob_big_diff

