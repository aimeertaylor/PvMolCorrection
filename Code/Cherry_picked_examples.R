################################################################################
# Cherry-picked examples of pretty data plus
# VHX_532 (an example of half sibs with a high posterior relapse probability)
# BPD_45 (a rare example of a possible reinfection with data on nine markers)
# All of which besides VHX_225_4 and VHX_532_4 (half-sibs) have common sense estimates.
################################################################################
rm(list = ls())
library(Pv3Rs)
cherries <- c("VHX_225", "VHX_457", "VHX_475", #"VHX_554", #"VHX_551", 
              "VHX_622") #, "VHX_532", "BPD_45" "VHX_541", "VHX_650")

load("../RData/results_Pv3Rs.RData")
marg <- sapply(cherries, function(pid) {
  if(is.na(Uniform_joint[pid])) {
    Uniform_pairwise[[pid]]
  } else {
    Uniform_joint[[pid]][["marg"]]
  }
})
state_names <- c("Recrud.", "Relapse", "Reinf.")
episodes <- unlist(sapply(marg, function(x) c(1, as.numeric(rownames(x)))))
probs <- unlist(sapply(marg, function(x) apply(x, 1, max))) # Get max marginal estimates 
states <- unlist(sapply(marg, function(x) apply(x, 1, function(z) state_names [which.max(z)])))
main_mar <- c(5,2.6,1,2.6) # Margin around main plot
plot_data(ys = ys_VHX_BPD[cherries], fs = fs_VHX_BPD, mar = main_mar, marker.annotate = F)
par(fig = c(0,1,0,1), mar = main_mar) # Reset before text annotation (important)
#box() # Toggle on and off to see plot limits for text() — not enacted in margin
text(y = rep(0.21, length(probs)), adj = 1, # right justify 
     x = seq(0.009, 1-0.009, length.out = length(episodes))[-which(episodes == 1)], 
     labels = paste0(states, " ", 100*round(probs,2), "%"), 
     cex = 0.7, srt = 90)

# For website: 
png("../Figures/cherry_picked.png")
plot_data(ys = ys_VHX_BPD[cherries], fs = fs_VHX_BPD, marker.annotate = F)
dev.off()


