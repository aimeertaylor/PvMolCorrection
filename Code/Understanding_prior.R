################################################################################
# When writing Compute_results_new.R, I noticed that the prior estimates plugged
# into the genetic model (i.e., the time-to-event posterior means) do not sum to
# one and, based on re-inspection of post_prob_CLI, these un-normalised
# estimates were used to compute the results in Taylor & Watson 2019 because
# post_prob_CLI only checks summation to one before taking logs if p_pop_ind =
# TRUE. As such, I saved the un-normalised and normalised prior estimates for
# further interrogation here. However, the differences are very small (for
# reinfection, at most 0.00015 in absolute terms and 0.2% in relative terms) so
# they likely have little effect.
################################################################################
rm(list = ls())

# Get episode IDs for episodes with genetic-based estimates
load("../../RData/marg_results_Pv3Rs.RData")
epIDs <- rownames(TimeToEvent_Pv3Rs)

# Load and trim priors to only episodes with genetic-based estimates
load("../../RData/prior_estimates.RData")
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



