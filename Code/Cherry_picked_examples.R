################################################################################
# Cherry-picked examples of pretty data plus
# VHX_532 (an example of half sibs with a high posterior relapse probability)
# BPD_45 (a rare example of a possible reinfection with data on nine markers)
# All of which besides VHX_225_4 and VHX_532_4 (half-sibs) have common sense estimates.
################################################################################
rm(list = ls())
library(Pv3Rs)
Fig <- TRUE
cherries <- c("VHX_225", "VHX_475", #"VHX_554", "VHX_551",
              "VHX_622", "VHX_532") #"BPD_45", "VHX_541", "VHX_650")

load("../RData/results_Pv3Rs.RData")
marg <- sapply(cherries, function(pid) {
  if(is.na(Uniform_joint[pid])) {
    Uniform_pairwise[[pid]]
  } else {
    Uniform_joint[[pid]][["marg"]]
  }
})

state_names <- c("C", "L", "I")
episodes <- unlist(sapply(cherries, function(pid) paste(pid, names(ys_VHX_BPD[[pid]]), sep = "_")))
probs <- unlist(sapply(marg, function(x) apply(x, 1, max))) # Get max marginal estimates 
states <- unlist(sapply(marg, function(x) apply(x, 1, function(z) state_names[which.max(z)])))
main_mar <- c(1.5, 3.5, 1.5, 3.5) # Margin around main plot
text <- rep("", length(episodes)) 
names(text) <- episodes
text[gsub(".", "_", names(probs), fixed = T)] <- paste0(states, " ", 100*round(probs,2))
z <- 1/(2*length(episodes)) # See text x placement
no_episodes_per_pid <- sapply(ys_VHX_BPD[cherries], length)
text_col <- unlist(sapply(1:length(cherries), function(i){
  if(i %% 2 == 1) {
    rep("white", no_episodes_per_pid[i])
  } else {
    rep("black", no_episodes_per_pid[i])
  }}))


if(Fig) png("../Figures/cherry_picked.png", res = 300, width = 10, height = 7, units = "in")
par(fig = c(0,1,0.2+0.01,1), mar = main_mar) # Reset before text annotation (important)
plot_data(ys = ys_VHX_BPD[cherries], fs = fs_VHX_BPD, mar = main_mar, marker.annotate = F)
par(fig = c(0,1,0.2+0.01,1), mar = main_mar) # Reset before text annotation (important)
text(y = rep(0.01, length(episodes)), x = seq(z, 1-z, length.out = length(episodes)), 
     labels = text, cex = 0.6, col = text_col)
if(Fig) dev.off()



