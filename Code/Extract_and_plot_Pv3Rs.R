# Set up and load new results
rm(list = ls())
library(Pv3Rs) # For simplex plots
load("../RData/results_Pv3Rs.RData") 
source("plot_VHXBPD_simplex.R") # wrapper to reduced copy-paste across scripts

#===============================================================================
# Uniform
#===============================================================================
# Create a list of all estimates and check no are NA
ls_Uniform_Pv3Rs <- sapply(Uniform_joint, function(x) x[[1]]) # Extract marginal 
ls_Uniform_Pv3Rs[is.na(Uniform_joint)] <- Uniform_pairwise[is.na(Uniform_joint)]
any(is.na(ls_Uniform_Pv3Rs))

# Create a named matrix from list
Uniform_Pv3Rs <- do.call(rbind, unname(ls_Uniform_Pv3Rs))
rownames(Uniform_Pv3Rs) <- unlist(sapply(names(ls_Uniform_Pv3Rs), function(pid) {
  paste(pid, rownames(ls_Uniform_Pv3Rs[[pid]]), sep = "_")
}))

# Extract joint / pairwise inference indicator and add to matrix
joint = unlist(sapply(names(ls_Uniform_Pv3Rs), function(pid) {
  nrec <- nrow(ls_Uniform_Pv3Rs[[pid]])
  rep(!is.na(Uniform_joint[pid]), nrec)}))
Uniform_Pv3Rs <- data.frame(Uniform_Pv3Rs, joint = joint)


#===============================================================================
# TimeToEvent
#===============================================================================
# Create a list of all estimates and check no are NA
ls_TimeToEvent_Pv3Rs <- sapply(TimeToEvent_joint, function(x) x[[1]]) # Extract marginal 
ls_TimeToEvent_Pv3Rs[is.na(TimeToEvent_joint)] <- TimeToEvent_pairwise[is.na(TimeToEvent_joint)]
any(is.na(ls_TimeToEvent_Pv3Rs))

# Create a named matrix from list
TimeToEvent_Pv3Rs <- do.call(rbind, unname(ls_TimeToEvent_Pv3Rs))
rownames(TimeToEvent_Pv3Rs) <- unlist(sapply(names(ls_TimeToEvent_Pv3Rs), function(pid) {
  paste(pid, rownames(ls_TimeToEvent_Pv3Rs[[pid]]), sep = "_")
}))

# Extract joint / pairwise inference indicator and add to matrix
joint = unlist(sapply(names(ls_TimeToEvent_Pv3Rs), function(pid) {
  nrec <- nrow(ls_TimeToEvent_Pv3Rs[[pid]])
  rep(!is.na(TimeToEvent_joint[pid]), nrec)}))
TimeToEvent_Pv3Rs <- data.frame(TimeToEvent_Pv3Rs, joint = joint)


#===============================================================================
# Save marginal results
#===============================================================================
save(Uniform_Pv3Rs, TimeToEvent_Pv3Rs, file = "../RData/marg_results_Pv3Rs.RData")

#===============================================================================
# Plot estimates by PMQ received or not
#===============================================================================
png("../Figures/simplex_Pv3Rs.png", width = 9, height = 9, units = "in", res = 300)
Uniform_xy <- apply(Uniform_Pv3Rs, 1, function(x) project2D(x[1:3]))
TimeToEvent_xy <- apply(TimeToEvent_Pv3Rs, 1, function(x) project2D(x[1:3]))
Uniform_xy <- rbind(Uniform_xy, joint = Uniform_Pv3Rs[,"joint"])
TimeToEvent_xy <- rbind(TimeToEvent_xy, joint = TimeToEvent_Pv3Rs[,"joint"])
plot_VHXBPD_simplex(Uniform_xy, TimeToEvent_xy)
dev.off()

# Extract data and estimates for recurrences whose recurrent state probabilities 
# could not be estimated using the prototype
unestimatable <- c( "VHX_239_2","VHX_461_2","VHX_39_2","VHX_52_2","VHX_583_2","VHX_33_2")
Uniform_Pv3Rs[unestimatable,]
TimeToEvent_Pv3Rs[unestimatable,]

# Check the other episodes of VHX_33 and VHX_583:
Uniform_Pv3Rs[grepl("VHX_33_", rownames(Uniform_Pv3Rs)),] # All to be expected. 
Uniform_Pv3Rs[grepl("VHX_583_", rownames(Uniform_Pv3Rs)),] # As expected


