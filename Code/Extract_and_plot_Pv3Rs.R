# Wrangle and plot new marginal results

# Set up and load new results
rm(list = ls())
load("../RData/results_Pv3Rs.RData") 
par(mfrow = c(2,2))

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

# Extract joint / pairwise indicator and add to matrix
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

# Extract joint / pairwise indicator and add to matrix
joint = unlist(sapply(names(ls_TimeToEvent_Pv3Rs), function(pid) {
  nrec <- nrow(ls_TimeToEvent_Pv3Rs[[pid]])
  rep(!is.na(TimeToEvent_joint[pid]), nrec)}))
TimeToEvent_Pv3Rs <- data.frame(TimeToEvent_Pv3Rs, joint = joint)

save(Uniform_Pv3Rs, TimeToEvent_Pv3Rs, file = "../RData/marg_results_Pv3Rs.RData")

#===============================================================================
# Plot estimates by study 
#===============================================================================
# Extract study indicator and check
BPD <- grepl("BPD", rownames(Uniform_Pv3Rs))
all(BPD == grepl("BPD", rownames(TimeToEvent_Pv3Rs)))

# Project
Uniform_xy <- apply(Uniform_Pv3Rs, 1, function(x) project2D(x[1:3]))
TimeToEvent_xy <- apply(TimeToEvent_Pv3Rs, 1, function(x) project2D(x[1:3]))
Uniform_LI_log <- Uniform_xy["y",] < -0.27 
TimeToEvent_LI_log <- Uniform_xy["y",] < -0.27 

# Plot BPD Uniform
plot_simplex(); title(main = "BPD Uniform", line = -2)
for(i in which(BPD)){
  points(x = Uniform_xy["x",i], y = Uniform_xy["y",i], 
         col = adjustcolor("black", alpha.f = 0.3), 
         pch = 17-Uniform_Pv3Rs[i,"joint"])}
lines(density(Uniform_xy["x",which(BPD)][Uniform_LI_log[which(BPD)]], 
              bw = "SJ", adjust = 5))

# Plot BPD TimeToEvent
plot_simplex(); title(main = "BPD Time-to-event", line = -2)
for(i in which(BPD)){
  points(x = TimeToEvent_xy["x",i], y = TimeToEvent_xy["y",i], 
         col = adjustcolor("black", alpha.f = 0.3), 
         pch = 17-TimeToEvent_Pv3Rs[i,"joint"])}
lines(density(TimeToEvent_xy["x",which(BPD)][TimeToEvent_LI_log[which(BPD)]], 
              bw = "SJ", adjust = 8))


# Plot VHX Uniform
plot_simplex(); title(main = "VHX Uniform", line = -2)
for(i in which(!BPD)){
  points(x = Uniform_xy["x",i], y = Uniform_xy["y",i], 
         col = adjustcolor("black", alpha.f = 0.3), 
         pch = 17-Uniform_Pv3Rs[i,"joint"])}
lines(density(Uniform_xy["x",which(!BPD)][Uniform_LI_log[which(!BPD)]], 
              bw = "SJ", adjust = 5))
legend(x = 0.1, y = 0.5, bty = "n", pch = 20, col = 9:8, inset = 0.01,
       legend = c("jointly modelled data","pairwise modelled data"))

# Plot VHX TimeToEvent
plot_simplex(); title(main = "VHX Time-to-event", line = -2)
for(i in which(!BPD)){
  points(x = TimeToEvent_xy["x",i], y = TimeToEvent_xy["y",i], 
         col = adjustcolor("black", alpha.f = 0.3), 
         pch = 17-TimeToEvent_Pv3Rs[i,"joint"])}
lines(density(TimeToEvent_xy["x",which(!BPD)][TimeToEvent_LI_log[which(!BPD)]], 
              bw = "SJ", adjust = 200))
legend(x = 0.1, y = 0.5, bty = "n", pch = 20, col = 9:8, inset = 0.01,
       legend = c("jointly modelled data","pairwise modelled data"))

par(mfrow = c(2,2))
plot(density(Uniform_xy["x",which(BPD)][Uniform_LI_log[which(BPD)]], 
              bw = "SJ", adjust = 5))
plot(density(TimeToEvent_xy["x",which(BPD)][TimeToEvent_LI_log[which(BPD)]], 
              bw = "SJ", adjust = 8))
plot(density(Uniform_xy["x",which(!BPD)][Uniform_LI_log[which(!BPD)]], 
              bw = "SJ", adjust = 5))
plot(density(TimeToEvent_xy["x",which(!BPD)][TimeToEvent_LI_log[which(!BPD)]], 
              bw = "SJ", adjust = 200))
