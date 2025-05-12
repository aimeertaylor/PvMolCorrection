################################################################################
# When writing Compute_results_new.R, I noticed that the prior estimates plugged
# into the genetic model (i.e., the time-to-event posterior means) do not sum to
# one and, based on re-inspection of post_prob_CLI, these unnormalised estimates
# were used to compute the old results, because post_prob_CLI only checks
# summation to one before taking logs if p_pop_ind = TRUE. As such, I saved the
# unnormalised and normalised prior estimates for further interrogation here.
#
# In summary, the differences are small (for reinfection, at most 0.00015 in
# absolute terms and 0.2% in relative terms) so they likely have little effect.
################################################################################
rm(list = ls())
library(Pv3Rs)

# Load prior estimates 
load("~/Dropbox/Vivax_VHXBPD_reanalysis/RData/prior_estimates.RData")

# Load results_new to get the episode IDs for genetic-based posterior estimates
load("~/Dropbox/Vivax_VHXBPD_reanalysis/RData/results_new.RData")
epIDs <- rownames(do.call(rbind, sapply(results_TimeToEvent, function(x) x["marg"])))

# Trim priors to only episodes with genetic-based posterior estimates
prior <- prior[epIDs, ]
prior_unnorm <- prior_unnorm[epIDs, ]

# Find discrepant priors
diff_ind <- !apply(prior == prior_unnorm, 1, all)

# Most discrepant 
diff_big <- epIDs[which.max(abs(prior[,"I"] - prior_unnorm[,"I"]))]

# Prior discrepancy scatter plots: 
plot(x = prior[diff_ind,"I"], 
     y = prior_unnorm[diff_ind,"I"], 
     pch = 4, cex = 0.5) # Not evident on 0-1 scale
points(x = prior[diff_big,"I"], 
       y = prior_unnorm[diff_big,"I"], 
       pch = 1, cex = 2, col = "red") 

plot(x = prior[diff_ind,"I"], 
     y = prior[diff_ind,"I"]-prior_unnorm[diff_ind,"I"],
     pch = 4, cex = 0.5) # Diff with P(I)
points(x = prior[diff_big,"I"], 
       y = prior[diff_big,"I"]-prior_unnorm[diff_big,"I"],
       pch = 1, cex = 2, col = "red") # Diff with P(I)

plot(x = prior[diff_ind,"I"], 
     y = (prior[diff_ind,"I"]-prior_unnorm[diff_ind,"I"])/prior[diff_ind,"I"],
     pch = 4, cex = 0.5) # Scaled diff with P(I)
points(x = prior[diff_big,"I"], 
     y = (prior[diff_big,"I"]-prior_unnorm[diff_big,"I"])/prior[diff_big,"I"],
     pch = 1, cex = 2, col = "red") # Scaled diff with P(I)



