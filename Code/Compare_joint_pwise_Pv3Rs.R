################################################################################
# Joint vs pairwise probabilities
# For loop suboptimally coded
# Something wrong with 419: compute_posterior(ys_VHX_BPD[["VHX_419"]] , fs_VHX_BPD, progress.bar = FALSE)
# recovers certain relapse for episode 5 and low relapse for episode 6!
################################################################################
rm(list = ls())
library(Pv3Rs) # For plot_data
load("../RData/results_Pv3Rs.RData")
big_diff <- 0.25 # probability difference considered large
states <- c(Recrudescence = "C", Reinfection = "I", Relapse = "L")
priors <- c("Uniform", "TimeToEvent")
pids_big_diff <- list()
eids_big_diff <- list()
Figs <- TRUE
marg_pwise <- list() # Store for extracted marginal probabilities
marg_joint <- list() # Store for extracted marginal probabilities

if(Figs) png("../Figures/compare_joint_vs_pwise.png", 
             width = 7, height = 10, units = "in", res = 300)
par(mfcol = c(3,2))

for(prior in priors) {
  
  # Load estimates
  if(prior == "Uniform"){
    joint <- Uniform_joint
    pairwise <- Uniform_pairwise
  } else {
    joint <- TimeToEvent_joint
    pairwise <- TimeToEvent_pairwise
  }
  
  # Extract marginal probabilities
  joint <- sapply(joint, function(x) x[[1]])
  
  # Tag pids that were not analysed jointly
  log_NA <- sapply(joint, function(x) any(is.na(x[[1]])))
  sum(log_NA) # Not modelled jointly 
  
  # Remove pids that were not analysed jointly
  joint <- joint[!log_NA]
  pwise <- pairwise[!log_NA]
  
  # Create matrices of estimates and summary features 
  unlst_joint <- do.call(rbind, joint)
  unlst_pwise <- do.call(rbind, pwise)
  rownames(unlst_joint) <- unlist(sapply(names(joint), function(pid) paste(pid, rownames(joint[[pid]]), sep = "_")))
  rownames(unlst_pwise) <- unlist(sapply(names(pwise), function(pid) paste(pid, rownames(pwise[[pid]]), sep = "_")))
  n_rec <- sapply(pwise, nrow) # Number of typed recurrences per pid
  unlst_n_rec <- rep(n_rec, n_rec) # Number of typed recurrences per pid per rec

  # Check people with only one recurrence have more-or-less identical estimates
  check <- max(unlst_joint[unlst_n_rec == 1,] - 
                 unlst_pwise[unlst_n_rec == 1,]) < .Machine$double.eps^0.5
  if(!check) stop("Discrepant estimates despite only one recurrence")
  
  # Plot estimates 
  for(s in states){
    
    plot(x = unlst_joint[unlst_n_rec > 1, s], 
         y = unlst_pwise[unlst_n_rec > 1, s], 
         xlim = c(0,1), ylim = c(0,1),
         ylab = "Pairwise", xlab = "Joint", pch = 20,
         bty = "n", main = sprintf("%s: %s prior", names(which(states == s)), prior))
    abline(a = 0, b = 1, lty = "dotted")
    
    # Extract and annotate big differences per recurrent state
    diffs <- abs(unlst_joint[unlst_n_rec > 1, s] - unlst_pwise[unlst_n_rec > 1, s]) 
    big_diffs <- names(diffs[which(diffs > big_diff)])
    if(length(big_diffs) > 0) {
      text(y = unlst_pwise[unlst_n_rec > 1, s][big_diffs], 
           x = unlst_joint[unlst_n_rec > 1, s][big_diffs], 
           labels = big_diffs, cex = 0.5,
           pos = if(s == "I") c(3,1,3,4,4,2) else c(1,3,1,2,2,4))
    }
  } 
  
  # Get differences per recurrencee
  big_diff_log <- apply(unlst_joint[unlst_n_rec > 1, ] - 
                          unlst_pwise[unlst_n_rec > 1, ], 1, 
                        function(x) any(abs(x) > big_diff))
  
  # Extract pids and eids for discrepant estimates
  pids_big_diff[[prior]] <- rep(names(n_rec[n_rec > 1]), n_rec[n_rec > 1])[big_diff_log]
  eids_big_diff[[prior]] <- names(big_diff_log)[big_diff_log]
  
  writeLines(sprintf("%s participants with data modelled jointly: %s recurrences
%s participants with data modelled jointly and more than one recurrence: %s recurrences
%s recurrences in %s people differ by more than %s percent probability", 
                     sum(!log_NA), nrow(unlst_joint), 
                     sum(n_rec > 1), sum(unlst_n_rec > 1),
                     sum(big_diff_log), length(pids_big_diff[[prior]]),
                     big_diff))
  
  # Discrepant estimates
  print(unlst_joint[unique(unlist(eids_big_diff[prior])),, drop = F])
  print(unlst_pwise[unique(unlist(eids_big_diff[prior])),, drop = F])
  
  # Store
  marg_pwise[[prior]] <- unlst_pwise
  marg_joint[[prior]] <- unlst_joint
  
}
if(Figs) dev.off()

#===============================================================================
# MS Figure
#===============================================================================
pids <- unique(unlist(pids_big_diff)) # All pids with a probability discrepant recurrence
unlist(eids_big_diff) == unique(unlist(eids_big_diff)) # No duplicate pids 

# Extract all episodes for pids with big difference
episodes <- unname(unlist(sapply(pids, function(pid) {
  sapply(names(ys_VHX_BPD[[pid]]), function(epi) {
    paste(pid, epi, sep = "_")})})))

# Extract probabilities for all recurrences in pids with at a discrepant recurrence
all_recs_pids_big_diff_Uniform <- unlist(sapply(pids_big_diff$Uniform, function(pid) paste(pid, names(ys_VHX_BPD[[pid]])[-1], sep = "_")))
all_recs_pids_big_diff_TimeToEvent <- unlist(sapply(pids_big_diff$TimeToEvent, function(pid) paste(pid, names(ys_VHX_BPD[[pid]])[-1], sep = "_")))

LorC_pwise <- 1-(rbind(marg_pwise[["Uniform"]][all_recs_pids_big_diff_Uniform,], 
                marg_pwise[["TimeToEvent"]][all_recs_pids_big_diff_TimeToEvent,,drop=F])[,"I"])
LorC_joint <- 1-(rbind(marg_joint[["Uniform"]][all_recs_pids_big_diff_Uniform,],
                 marg_joint[["TimeToEvent"]][all_recs_pids_big_diff_TimeToEvent,,drop=F])[,"I"])

no_episodes_per_pid <- sapply(pids, function(pid) length(ys_VHX_BPD[[pid]]))
text <- text_pwise <- text_joint <- rep("", length(episodes)) 
names(text_pwise) <- names(text_joint) <- episodes
text_pwise[names(LorC_pwise)] <- round(LorC_pwise, 2)*100
text_joint[names(LorC_joint)] <- round(LorC_joint, 2)*100
text_col <- unlist(sapply(1:length(pids), function(i){
  if(i %% 2 == 1) {
    rep("white", no_episodes_per_pid[i])
  } else {
    rep("black", no_episodes_per_pid[i])
  }}))
main_mar <- c(1.5, 3.5, 1.5, 3.5)
z <- 1/(2*length(episodes)) # See text x placement

# Plot for ms
if(Figs) png("../Figures/data_joint_vs_pwise.png", width = 10, height = 7, units = "in", res = 300) 
par(fig = c(0,1,0.2+0.01,1), mar = main_mar) # Important to call before and after plot_data
plot_data(ys = ys_VHX_BPD[pids], fs = fs_VHX_BPD, marker.annotate = F, mar = main_mar)
par(fig = c(0,1,0.2+0.01,1), mar = main_mar) # Important to call before and after plot_data
text(y = rep(0.04, length(episodes)),
     x = seq(z, 1-z, length.out = length(episodes)), 
     labels = text_pwise, cex = 0.6, col = "black")
text(y = rep(0.01, length(episodes)),
     x = seq(z, 1-z, length.out = length(episodes)), 
     labels = text_joint, cex = 0.6, col = text_col)
points(y = rep(-0.01, length(unlist(eids_big_diff))),
       x = seq(z, 1-z, length.out = length(episodes))[episodes %in% unlist(eids_big_diff)], 
       pch = 17, cex = 0.5) 
if(Figs) dev.off()

load("../RData/prior_estimates.RData")
prior["VHX_419_6",]
