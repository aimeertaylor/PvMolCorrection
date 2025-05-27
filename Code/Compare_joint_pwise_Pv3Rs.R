################################################################################
# Joint vs pairwise probabilities
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
  n_rec <- sapply(pwise, nrow) # Number of recurrences per pid
  unlst_n_rec <- rep(n_rec, n_rec) # Number of recurrences per pid per rec
  rownames(unlst_joint) <- rownames(unlst_pwise) <- # Add estimate IDs
    unlist(sapply(names(n_rec), function(x) paste(x, 1:n_rec[x]+1, sep = "_")))
  
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
                     sum(big_diff_log), length(pids_big_diff),
                     big_diff))
  
  # Discrepant estimates
  print(unlst_joint[unique(unlist(eids_big_diff[prior])),, drop = F])
  print(unlst_pwise[unique(unlist(eids_big_diff[prior])),, drop = F])
}
if(Figs) dev.off()

# Data plot 
if(Figs) png("../Figures/data_joint_vs_pwise.png", 
             width = 7, height = 7, units = "in", res = 300) 
plot_data(ys = ys_VHX_BPD[unique(unlist(pids_big_diff))], fs = fs_VHX_BPD)
if(Figs) dev.off()



