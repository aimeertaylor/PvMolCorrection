################################################################################
# Cherry-picked examples of pretty data plus
# VHX_532 (an example of half sibs with a high posterior relapse probability)
# BPD_45 (a rare example of a possible reinfection with data on nine markers)
# All of which besides VHX_225_4 and VHX_532_4 (half-sibs) have common sense estimates.
################################################################################
rm(list = ls())
library(Pv3Rs)
cherries <- c("VHX_225", "VHX_457", "VHX_475", "VHX_554", "VHX_551", 
              "VHX_622", "VHX_532", "BPD_45")# "VHX_541", "VHX_650")

load("../RData/results_Pv3Rs.RData")
marg <- sapply(cherries, function(pid) {
  if(is.na(Uniform_joint[pid])) {
    Uniform_pairwise[[pid]]
  } else {
    Uniform_joint[[pid]][["marg"]]
  }
})
state_names <- c("Recrudescence", "Relapse", "Reinfection")
episodes <- unlist(sapply(marg, function(x) c(1, as.numeric(rownames(x)))))
probs <- unlist(sapply(marg, function(x) apply(x, 1, max))) # Get max marginal estimates 
states <- unlist(sapply(marg, function(x) apply(x, 1, function(z) state_names [which.max(z)])))
plot_data(ys = ys_VHX_BPD[cherries], fs = fs_VHX_BPD)
text(x = rep(0.07, length(probs)), 
     y = seq(0.05, 0.95, length.out = length(episodes))[-which(episodes == 1)], 
     labels = paste0(states, " ", 100*round(probs,2), "%"), 
     cex = 0.7)

# For website: 
png("../Figures/cherry_picked.png")
plot_data(ys = ys_VHX_BPD[cherries], fs = fs_VHX_BPD, marker_annotate = F)
dev.off()



