################################################################################
# Script to explore inference as a function of markers Assumes pairwise
# inference: written with Emily's False Negatives in mind.
#
# Examples without suspected misspecification: 
# "VHX_52", "BPD_562", "VHX_239", "VHX_461"
# Examples with suspected misspecification:
# "VHX_56", "VHX_39", "BPD_577"
#
# Emily will need to replace:
# ys_VHX_BPD with ys from OPRA, IMPROV
# fs_VHX_BPD with fs from OPRA, IMPROV
#
# In example below we see that when inference is based on data that include
# PV.ms8, relapse probability generally drops to its lowest possible value
################################################################################
library(Pv3Rs)
rm(list = ls())
oldpar <- par(no.readonly = T) # To restore plotting parameters later
pid <- "VHX_56" # Participant id of interest
plot_data(ys_VHX_BPD[pid], fs = fs_VHX_BPD) # Plot Data
markers <- names(fs_VHX_BPD) # Names of markers
all_permutations <- gtools::permutations(repeats.allowed = F, # All possible permutations
                                            n = length(markers), 
                                            r = length(markers))
set.seed(1) # marker reproducible
n_permutations <- 10 # Number of permutations to explore must be <= nrow(marker_permutations)
permutations_to_plot <- sample(x = nrow(all_permutations), size = n_permutations, replace = F)

for(permutation in permutations_to_plot){
  
  marker_permute <- markers[all_permutations[permutation,]]
   
  results <- lapply(1:length(markers), function(m){ # for marker subsets
    
    # extract data on maker subset
    y <- lapply(ys_VHX_BPD[[pid]], function(y_epi) y_epi[marker_permute][1:m]) 
    # compute posterior probabilities
    z <- compute_posterior(y, fs = fs_VHX_BPD, return.RG = TRUE, return.logp = TRUE)
    # extract log likelihoods of relationship graphs RG
    llikes <- sapply(z$RGs, function(RG) RG$logp)
    # extract maximum likelihood (ML) RG 
    RG <- z$RGs[[which(abs(llikes - max(llikes)) < .Machine$double.eps^0.5)[1]]]
    # return probabilities and ML RGs 
    list(probs = z$marg, RG = RG_to_igraph(RG, determine_MOIs(y)))
  })
  
  posterior_probs <- sapply(results, function(x) x[["probs"]]) # extract probs
  
  # Plot evolution of reinfection probability
  graphics::layout(mat = matrix(c(rep(1, 9), 2:10), nrow = 2, byrow = T))
  par(mar = c(4,5,1,2), pty = "m")
  plot(y = posterior_probs[2,], x = 1:length(markers), 
       xlab = "Number of markers", ylab = "Probability of relapse", 
       ylim = c(0,1), bty = "n", type = "b", pch = 20)
  
  # Add ML RGs
  par(mar = rep(1.5,4), pty = "s")
  for(i in 1:length(results)){
    plot_RG(results[[i]][["RG"]], vertex.size = 20, vertex.label = NA)
    title(main = marker_permute[i]) # Number of markers
  }
  
  par(oldpar)
}
