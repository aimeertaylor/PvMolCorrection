################################################################################
# In this script old and new posterior estimates are compared
#
# Run this script after Compute_recurrent_states_new.R and before
# Compute_PQ_failure_rates_new.R --- it outputs Results_BPD needed for
# Compute_PQ_failure_rates_new.R 
#
# Within the old results, we can either compare to the median or the mean (they
# are not equal - see plots below). Besides supplementary figure 9, almost all
# figures and computations in Pooled_Analysis.Rmd and thus the Taylor&Watson
# 2019 are based on the median. This includes the PQ failure rate computations.
# As such, the old median is needed to internally check failure rate
# computation. The old mean is needed to compare with the new mean - the new
# median is computationally expensive to re-estimate, with a small return on
# investment, especially if new estimates are not included in a manuscript. To
# extract both median and mean estimates for Compute_PQ_failure_rates_new.R, I
# thus loop over median and mean. This is suboptimal coding --- would have been
# better to save both median and mean in Results_BPD --- but it does the job
# (the time investment to optimise is non-negligible for a small return).
# 
# To ensure the only difference between new and old results is the updated
# model, we use the old allele frequencies and compare new results with those
# that were directly computed in Taylor & Watson et al. 2019. Those that were
# directly computed in Taylor & Watson et al. 2019 exclude all those with data
# on more than three episodes (default argument Max_Eps = 3 of function
# post_prob_CLI). For three of nine patients whose data were analysed when
# EXCLUDE_COMPLEX = F (the option to exclude these patients was set to aid
# illustrative runs of the code), NAs were returned. That said, results
# generated using the time-to-event posterior estimates as prior estimates also
# differ because prior estimates plugged into the old model were un-normalised
# but are normalised here.
#
# To-do: 
# Make paths run on the internet
# Save illustrative example in another file
################################################################################
rm(list = ls())
library(Pv3Rs)
pardefault <- par()

# https://www.statology.org/euclidean-distance-in-r/
euclidean <- function(a, b) sqrt(sum((a - b)^2))

# Function to extract states whose posterior probability exceeds p
get_post_states <- function(x, p = 0.1) {
  paste(sort(colnames(new_Uniform)[which(x > p)]), collapse = "_")
}

# Load new results
load("~/Dropbox/Vivax_VHXBPD_reanalysis/RData/results_new.RData") 

# Extract marginal probabilities
new_Uniform <- do.call(rbind, sapply(results_Uniform, function(x) x["marg"]))
new_TimeToEvent <- do.call(rbind, sapply(results_TimeToEvent, function(x) x["marg"]))

# Load directly-estimated old results 
path <- "~/Documents/RecurrentVivax" # Path to old estimates
load(sprintf('%s/RData/GeneticModel/Including_Complex_Cases_Full_Posterior_Model_samples_Tagnostic.RData', path))
load(sprintf('%s/RData/GeneticModel/Including_Complex_Cases_Full_Posterior_Model_samples.RData', path))

# ==============================================================================
# Extract old median and mean results 
# ==============================================================================
par(mfrow = c(2,1))

# Small variation between mean and median using the uniform prior suggests
# variation due to time-to-event draws (which we can exactly reproduce)
# dominates variation due to posterior allele-frequency draws (which we cannot
# exactly reproduce). As such, generation of estimates for all 100 posterior
# time-to-event is not rendered worthless by the fact that we cannot recreate
# draws from the allele frequency posterior exactly (we didn't use set.seed in
# Pooled_Analysis.Rmd). Nonetheless, it is computationally expensive
plot(y = thetas_9MS$I, x = thetas_9MS$`I50%`, 
     main = "Time-to-event prior", 
     ylab = "mean", xlab = "median")
plot(y = thetas_9MS_Tagnostic$I, x = thetas_9MS_Tagnostic$`I50%`,
     main = "Uniform prior",
     ylab = "mean", xlab = "median")

# ==============================================================================
# Run for both old median and old mean
# ==============================================================================
for (comparator in c("median", "mean")) { 
  
  if (comparator == "mean") {
    old_Uniform <- thetas_9MS_Tagnostic[, c("C", "L", "I")] 
    old_TimeToEvent <- thetas_9MS[, c("C", "L", "I")]   
  } else {
    old_Uniform <- thetas_9MS_Tagnostic[, c("C50%", "L50%", "I50%")]
    old_TimeToEvent <- thetas_9MS[, c("C50%", "L50%", "I50%")]  
  }
  colnames(old_Uniform) <- c("C","L","I")
  colnames(old_TimeToEvent) <- c("C","L","I")
  
  
  # ==============================================================================
  # Check the same episodes are directly-estimated within new and old
  # ==============================================================================
  setequal(row.names(new_TimeToEvent), row.names(new_Uniform))
  setequal(row.names(old_TimeToEvent), row.names(old_Uniform))
  
  # ==============================================================================
  # Check for NA results in the old results 
  # ==============================================================================
  old_Uniform_NAs <- which(apply(is.na(old_Uniform), 1, any))
  old_TimeToEvent_NAs <- which(apply(is.na(old_TimeToEvent), 1, any))
  setequal(old_Uniform_NAs, old_TimeToEvent_NAs)
  
  # ==============================================================================
  # Remove NA results from old
  # ==============================================================================
  old_Uniform <- old_Uniform[-old_Uniform_NAs, ]
  old_TimeToEvent <- old_TimeToEvent[-old_TimeToEvent_NAs, ]
  
  # ==============================================================================
  # Extract comparable 
  # ==============================================================================
  comparable <- intersect(rownames(new_Uniform), rownames(old_Uniform))  
  incomparable_new <- rownames(new_Uniform)[!rownames(new_Uniform) %in% comparable]
  incomparable_old <- rownames(old_Uniform)[!rownames(old_Uniform) %in% comparable]
  
  # ==============================================================================
  # List incomparable / comparable
  # ==============================================================================
  print(incomparable_new)
  print(incomparable_old)
  
  writeLines(sprintf("The new model is able to directly compute posterior 
                   estimates for %s additional recurrences",
                   length(incomparable_new)))
  
  
  # Are the BPD all comparable? Yes! 
  if(setequal(rownames(old_Uniform)[grepl("BPD", rownames(old_Uniform))], 
              rownames(new_Uniform)[grepl("BPD", rownames(new_Uniform))])) {
    print("All BPD are comparable")
  } else {
    print("Some BPD incomparable")
  }
  
  # List comparable episodes 
  comparable_BPD <- comparable[grepl("BPD", comparable)]
  comparable_VHX <- comparable[grepl("VHX", comparable)]
  
  
  #===============================================================================
  # Extract and save BPD time-to-event estimates for Compute_PQ_failure_rates_new.R
  #===============================================================================
  Results_BPD <- as.data.frame(cbind(old = old_TimeToEvent[comparable_BPD, "I"], 
                                     new = new_TimeToEvent[comparable_BPD, "I"]))
  save(Results_BPD, file = 
         sprintf("~/Dropbox/Vivax_VHXBPD_reanalysis/RData/Results_BPD_%s.RData", comparator))
  
  #===============================================================================
  # Posterior estimates put weight on either one or two states
  #===============================================================================
  # L, I, L&I (Time to event), 
  # L, I, L&I, C, L&C (Uniform) 
  unique(apply(new_Uniform, 1, get_post_states))
  unique(apply(old_Uniform, 1, get_post_states))
  unique(apply(new_TimeToEvent, 1, get_post_states))
  unique(apply(old_TimeToEvent, 1, get_post_states))
  
  
  #===============================================================================
  # Compute differences 
  #===============================================================================
  # Project probabilities onto 2D simplex coordinates 
  xy0_u <- apply(old_Uniform[comparable,], 1, project2D)
  xy1_u <- apply(new_Uniform[comparable,], 1, project2D)
  xy0_t <- apply(old_TimeToEvent[comparable,], 1, project2D)
  xy1_t <- apply(new_TimeToEvent[comparable,], 1, project2D)
  
  # Compute distance on 2D simplex
  diffs_Uniform <- sapply(1:ncol(xy0_u), function(i) euclidean(xy0_u[,i], xy1_u[,i]))
  diffs_TimeToEvent <- sapply(1:ncol(xy0_t), function(i) euclidean(xy0_t[,i], xy1_t[,i]))
  
  # Plot histogram of differences
  par(mfrow = c(1,2))
  hist(diffs_Uniform)
  hist(diffs_TimeToEvent)
  
  # Get the names of the episodes that diff most
  big_diff <- 0.2
  big_diff_u <- colnames(xy0_u)[diffs_Uniform > big_diff]
  big_diff_t <- colnames(xy0_u)[diffs_TimeToEvent > big_diff]
  
  # Most divergent in BPD and VHX respectively
  big_diff_u_BPD <- intersect(big_diff_u, comparable_BPD)
  big_diff_t_BPD <- intersect(big_diff_t, comparable_BPD)
  big_diff_u_VHX <- intersect(big_diff_u, comparable_VHX)
  big_diff_t_VHX <- intersect(big_diff_t, comparable_VHX)
  
  # Save for IllustrativeExamples.R
  big_diff_eid <- union(big_diff_u, big_diff_t)
  save(big_diff_eid, file = sprintf("~/Dropbox/Vivax_VHXBPD_reanalysis/RData/big_diff_eid_%s.RData", comparator))
  save(old_TimeToEvent, old_Uniform, 
       file = sprintf("~/Dropbox/Vivax_VHXBPD_reanalysis/RData/results_old_%s.RData", comparator))
  
  #===============================================================================
  # Scatter plots of relapse probabilities  
  # 2 of 8 large deviations are undesirable: VHX_91_2 and VHX_56_2 
  # Another good example: VHX_91_3, but not as impacted because VHX_91_3 not polyclonal
  # Write a vignette around VHX_91_2 and VHX_91_3 and VHX_56_2, 
  # VHX_91 seems to be a true half, illustrating problem of new model assuming all 
  # siblings are half siblings. VHX_56_2 illustrates problem of new model relapse 
  # being more brittle to genotyping errors where it was not before
  #===============================================================================
  
  
  png(sprintf("~/Dropbox/Vivax_VHXBPD_reanalysis/Plots/Compare_scatter_%s.png", comparator))
  par(mfrow = c(2,2))
  plot(NULL, xlim = c(0,1), ylim = c(0,1), bty = "n", 
       main = "VHX: uniform prior", 
       xlab = "Old relapse posterior", 
       ylab = "New relapse posterior")
  abline(a = 0, b = 1, col = "lightgray")
  points(x = old_Uniform[comparable_VHX, "L"], 
         y = new_Uniform[comparable_VHX, "L"], 
         pch = 4, col = "cornflowerblue")
  if (length(setdiff(big_diff_t_VHX, big_diff_u_VHX)) > 0 ){
    text(x = old_Uniform[setdiff(big_diff_t_VHX, big_diff_u_VHX), "L"], 
         y = new_Uniform[setdiff(big_diff_t_VHX, big_diff_u_VHX), "L"], 
         labels = setdiff(big_diff_t_VHX, big_diff_u_VHX), 
         pos = 3, cex = 0.45, xpd = NA)}
  if (length(big_diff_u_VHX) > 0 ){
    text(x = old_Uniform[big_diff_u_VHX, "L"], 
         y = new_Uniform[big_diff_u_VHX, "L"], 
         labels = big_diff_u_VHX, 
         pos = 3, cex = 0.75, xpd = NA)}
  
  plot(NULL, xlim = c(0,1), ylim = c(0,1), bty = "n", 
       main = "VHX: time-to-event prior",
       xlab = "Old relapse posterior", 
       ylab = "New relapse posterior")
  abline(a = 0, b = 1, col = "lightgray")
  points(x = old_TimeToEvent[comparable_VHX, "L"], 
         y = new_TimeToEvent[comparable_VHX, "L"], 
         pch = 4, col = "cornflowerblue")
  if (length(setdiff(big_diff_u_VHX, big_diff_t_VHX)) > 0 ){
    text(x = old_TimeToEvent[setdiff(big_diff_u_VHX, big_diff_t_VHX), "L"], 
         y = new_TimeToEvent[setdiff(big_diff_u_VHX, big_diff_t_VHX), "L"], 
         labels = setdiff(big_diff_u_VHX, big_diff_t_VHX), 
         pos = 3, cex = 0.45, xpd = NA)}
  if (length(big_diff_t_VHX) > 0 ){
    text(x = old_TimeToEvent[big_diff_t_VHX, "L"], 
         y = new_TimeToEvent[big_diff_t_VHX, "L"], 
         labels = big_diff_t_VHX, 
         pos = 3, cex = 0.75, xpd = NA)}
  
  
  plot(NULL, xlim = c(0,1), ylim = c(0,1), bty = "n", 
       main = "BPD: Uniform prior", 
       xlab = "Old relapse posterior", 
       ylab = "New relapse posterior")
  abline(a = 0, b = 1, col = "lightgray")
  points(x = old_Uniform[comparable_BPD, "L"], 
         y = new_Uniform[comparable_BPD, "L"], 
         pch = 4, col = "cornflowerblue")
  if (length(setdiff(big_diff_t_BPD,  big_diff_u_BPD)) > 0){
    text(x = old_Uniform[setdiff(big_diff_t_BPD,  big_diff_u_BPD), "L"], 
         y = new_Uniform[setdiff(big_diff_t_BPD,  big_diff_u_BPD), "L"], 
         labels = setdiff(big_diff_t_BPD,  big_diff_u_BPD), 
         pos = 3, cex = 0.45, xpd = NA)}
  if (length(big_diff_u_BPD) > 0 ){
    text(x = old_Uniform[big_diff_u_BPD, "L"], 
         y = new_Uniform[big_diff_u_BPD, "L"], 
         labels = big_diff_u_BPD, 
         pos = 3, cex = 0.75, xpd = NA)}
  
  plot(NULL, xlim = c(0,1), ylim = c(0,1), bty = "n",
       main = "Time-to-event prior",
       xlab = "Old relapse posterior", 
       ylab = "New relapse posterior")
  abline(a = 0, b = 1, col = "lightgray")
  points(x = old_TimeToEvent[comparable_BPD, "L"], 
         y = new_TimeToEvent[comparable_BPD, "L"], 
         pch = 4, col = "cornflowerblue")
  if (length(setdiff(big_diff_u_BPD, big_diff_t_BPD)) > 0 ){
    text(x = old_TimeToEvent[setdiff(big_diff_u_BPD, big_diff_t_BPD), "L"], 
         y = new_TimeToEvent[setdiff(big_diff_u_BPD, big_diff_t_BPD), "L"], 
         labels = setdiff(big_diff_u_BPD, big_diff_t_BPD), 
         pos = 3, cex = 0.45, xpd = NA)}
  if (length(big_diff_t_BPD) > 0){
    text(x = old_TimeToEvent[big_diff_t_BPD, "L"], 
         y = new_TimeToEvent[big_diff_t_BPD, "L"], 
         labels = big_diff_t_BPD, 
         pos = 3, cex = 0.75, xpd = NA)}
  
  dev.off()
  
  #===============================================================================
  # Simplex plots Highlight the fact that posterior probability falls on one or
  # two states However, because posterior probability falls on one or two states,
  # simplex plots are not very useful (it wouldn't necessarily be easier to see
  # vectors on L-I plane or on L-C plane - everything would fall on the
  # diagonal).
  #===============================================================================
  
  png(sprintf("~/Dropbox/Vivax_VHXBPD_reanalysis/Plots/Compare_simplex_%s.png", comparator))
  # Prior discrepancy simplex plot (not evident)
  par(mar = c(0,0,0,0), mfrow = c(2,2))
  prior_Uniform <- project2D(rep(1/3,3))
  
  plot_simplex(v_labels = colnames(old_Uniform))
  text(x = prior_Uniform["x"], y = prior_Uniform["y"], labels = "BPD Uniform prior", pos = 3)
  points(x = xy1_u["x",comparable_BPD], y = xy1_u["y",comparable_BPD], 
         col = "cornflowerblue", cex= 2, pch = 20)
  if (length(big_diff_u_BPD) > 0 ){
    arrows(x0 = xy0_u["x",big_diff_u], x1 = xy1_u["x",big_diff_u],
           y0 = xy0_u["y",big_diff_u], y1 = xy1_u["y",big_diff_u], 
           length = 0.05, col = "red")
  }
  
  plot_simplex(v_labels = colnames(old_TimeToEvent))
  text(x = prior_Uniform["x"], y = prior_Uniform["y"], labels = "BPD Time-to-event prior", pos = 3)
  points(x = xy1_t["x",comparable_BPD], y = xy1_t["y",comparable_BPD], 
         col = "cornflowerblue", cex= 2, pch = 20)
  if (length(big_diff_t_BPD) > 0 ){
    arrows(x0 = xy0_t["x",big_diff_t_BPD], x1 = xy1_t["x",big_diff_t_BPD],
           y0 = xy0_t["y",big_diff_t_BPD], y1 = xy1_t["y",big_diff_t_BPD], 
           length = 0.05, col = "red")
  }
  
  plot_simplex(v_labels = colnames(old_Uniform))
  text(x = prior_Uniform["x"], y = prior_Uniform["y"], labels = "VHX Uniform prior", pos = 3)
  points(x = xy1_u["x",comparable_VHX], y = xy1_u["y",comparable_VHX], 
         col = "cornflowerblue", cex= 2, pch = 20)
  if (length(big_diff_u_VHX) > 0 ){
    arrows(x0 = xy0_u["x",big_diff_u], x1 = xy1_u["x",big_diff_u],
           y0 = xy0_u["y",big_diff_u], y1 = xy1_u["y",big_diff_u], 
           length = 0.05, col = "red")
  }
  
  plot_simplex(v_labels = colnames(old_TimeToEvent))
  text(x = prior_Uniform["x"], y = prior_Uniform["y"], labels = "VHX Time-to-event prior", pos = 3)
  points(x = xy1_t["x",comparable_VHX], y = xy1_t["y",comparable_VHX], 
         col = "cornflowerblue", cex= 2, pch = 20)
  if (length(big_diff_t_VHX) > 0 ){
    arrows(x0 = xy0_t["x",big_diff_t_VHX], x1 = xy1_t["x",big_diff_t_VHX],
           y0 = xy0_t["y",big_diff_t_VHX], y1 = xy1_t["y",big_diff_t_VHX], 
           length = 0.05, col = "red")
  }
  dev.off()
}

# Restore plotting margins
par(mar = pardefault$mar)



