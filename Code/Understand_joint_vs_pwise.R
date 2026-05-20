################################################################################
# By looking at individual cases, this script digs deeper into the general
# observation that both reinfection and recrudescence have slightly elevated
# probabilities given joint vs approximate joint inf
################################################################################
rm(list = ls())
library(Pv3Rs)
source("compute_posterior_approxjoint.R")
load("../RData/marg_joint_pwise.RData") # From Compare_joint_pwise_Pv3Rs.R

#===============================================================================
# Synthetic data show that
#
# probability of per-recurrence reinfection increases with the episode count
# within large graphs.
#
# probability of per-recurrence recrudescence is constant across episodes within
# large graphs but that constant increases with graph size.
#
# approximate joint inference cannot capture either of these things (these
# results are not in Understanding posterior probabilities)
#===============================================================================
fs <- list(m1 = c('1' = 0.2, '2' = 0.2, '3' = 0.2, '4' = 0.2, '5' = 0.2),
           m2 = c('1' = 0.2, '2' = 0.2, '3' = 0.2, '4' = 0.2, '5' = 0.2),
           m3 = c('1' = 0.2, '2' = 0.2, '3' = 0.2, '4' = 0.2, '5' = 0.2))

# Synthetic reinfection data: 
yhet <- list("0" = list(m1 = "1", m2 = "1", m3 = "1"),
             "1" = list(m1 = "2", m2 = "2", m3 = "2"),
             "2" = list(m1 = "3", m2 = "3", m3 = "3"),
             "3" = list(m1 = "4", m2 = "4", m3 = "4"),
             "4" = list(m1 = "5", m2 = "5", m3 = "5"))
plot_data(list("synthetic reinfections" = yhet), fs)
suppressMessages(compute_posterior(yhet, fs, progress.bar = FALSE))$marg
suppressMessages(compute_posterior_approxjoint(yhet, fs))

# Synthetic recrudescence data: 
yhoms <- list(y1 = list("0" = list(m1 = "1", m2 = "1", m3 = "1"),
                        "1" = list(m1 = "1", m2 = "1", m3 = "1")), 
              y2 = list("0" = list(m1 = "1", m2 = "1", m3 = "1"),
                        "1" = list(m1 = "1", m2 = "1", m3 = "1"),
                        "2" = list(m1 = "1", m2 = "1", m3 = "1")),
              y3 = list("0" = list(m1 = "1", m2 = "1", m3 = "1"),
                        "1" = list(m1 = "1", m2 = "1", m3 = "1"),
                        "2" = list(m1 = "1", m2 = "1", m3 = "1"),
                        "3" = list(m1 = "1", m2 = "1", m3 = "1")),
              y4 = list("0" = list(m1 = "1", m2 = "1", m3 = "1"),
                        "1" = list(m1 = "1", m2 = "1", m3 = "1"),
                        "2" = list(m1 = "1", m2 = "1", m3 = "1"),
                        "3" = list(m1 = "1", m2 = "1", m3 = "1"),
                        "4" = list(m1 = "1", m2 = "1", m3 = "1")))
sapply(yhoms, function(yhom) suppressMessages(compute_posterior(yhom, fs))$marg)
sapply(yhoms, function(yhom) suppressMessages(compute_posterior_approxjoint(yhom, fs)))


#===============================================================================
# At a study level, VHX BPD show that, in general, insight based on synthetic data translates 
# to real data 
#===============================================================================

# Differences between joint and approximate joint given probable reinfection:
# Increases with episode count as predicted by synthetic data above
diffs <- marg_joint[["Uniform"]][, "I"] - marg_pwise[["Uniform"]][, "I"]
reinf_eids <- names(which(marg_joint[["Uniform"]][, "I"] > 0.5)) # Important to condition on probable reinfection
diffs_reinf <- diffs[reinf_eids]
eid_no_reinf <- as.numeric(do.call(rbind, strsplit(names(diffs_reinf), split = "_"))[,3])
plot(x = eid_no_reinf, y = diffs_reinf, ylim = c(-0.25, 0.25), # Exclude outlies with diff > 0.25
     xlab = "Episode number", ylab = "Absolute difference between joint and approximate joint",
     main = "Probable reinfection")

# Differences between joint and approximate joint given probable recrudescence:
diffs <- marg_joint[["Uniform"]][, "C"] - marg_pwise[["Uniform"]][, "C"]
recru_eids <- names(which(marg_joint[["Uniform"]][, "C"] > 0.5)) 
diffs_recru <- diffs[recru_eids]
pids <- apply(do.call(rbind, strsplit(names(diffs_recru), split = "_"))[,-3], 1, paste, collapse = "_")
maxepi_no_recru <- sapply(ys_VHX_BPD[pids], length)
plot(x = maxepi_no_recru, y = diffs_recru, ylim = c(-0.25, 0.25), # Exclude outlies with diff > 0.25
     xlab = "Maimum episode number", ylab = "Absolute difference between joint and approximate joint", 
     main = "Probable recrudescence")
abline(lm(diffs_recru ~ maxepi_no_recru)) # Linear increase

#===============================================================================
# At a per-participant level, VHX BPD show that, there are cases that conform to
# the insight gleamed from synthetic data — e.g., for refection see VHX_621 (some matches and a
# clear relapse) and BPD_455 (some matches) and BPD_234 (no matches); for recrudescence see BPD_347 — but also
# cases where reinfection appears elevated because the likelihood on the all
# sibling graph is zero (VHX_66)
#===============================================================================
# Choose pid and plot data:
pid <- "VHX_66"
plot_data(ys_VHX_BPD[pid], fs = sapply(fs_VHX_BPD, sort, decreasing = T), 
          marker.annotate = F)

# Compute recurrence state probabilities:
joint <- compute_posterior(ys_VHX_BPD[[pid]], fs = fs_VHX_BPD, return.logp = T)          
aprox <- compute_posterior_approxjoint(ys_VHX_BPD[[pid]], fs = fs_VHX_BPD)

# Inspect recurrence state probabilities:
joint$marg
aprox

if (pid == "VHX_66") {
  suppressMessages(compute_posterior(ys_VHX_BPD[[pid]][c(1,2)], fs = fs_VHX_BPD))$marg
  suppressMessages(compute_posterior(ys_VHX_BPD[[pid]][c(1,3)], fs = fs_VHX_BPD))$marg
  suppressMessages(compute_posterior(ys_VHX_BPD[[pid]][c(2,3)], fs = fs_VHX_BPD))$marg
  
  RG_liks <- sapply(joint$RGs, function(x) x$logp)
  RG_inds <- which(!is.infinite(RG_lik))
  mois <- determine_MOIs(ys_VHX_BPD[[pid]])
  par(mfrow = c(3,2))
  sapply((1:length(RG_liks)), function(ind) {
    RG <- joint$RGs[[ind]]
    plot_RG(RG_to_igraph(RG, mois), edge.curved = 0.5)
    title(RG_liks[ind])
    box()})
  
  # For VHX_66, I think if the likelihood of the all sibling graph were not -Inf,
  # the probabilities of IL and II would be equal whereas:
  sort(joint$joint, decreasing = T)
  abs(joint$marg[,"I"] - aprox[,"I"])  
}

