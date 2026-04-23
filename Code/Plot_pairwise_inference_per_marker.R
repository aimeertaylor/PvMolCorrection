################################################################################
# Script to explore inference as a function of markers. Assumes pairwise
# inference: written with Emily's false negatives in mind.
#
# Examples without suspected misspecification: 
# "VHX_52", "BPD_562", "VHX_239", "VHX_461"
# Examples with suspected misspecification:
# "VHX_56", "VHX_39", "BPD_577"
#
# Emily will need to replace:
# ys_VHX_BPD with ys from OPRA, IMPROV
# fs_VHX_BPD with fs from OPRA, IMPROV
################################################################################
library(Pv3Rs)
rm(list = ls())
set.seed(1) # marker reproducible
pid <- "VHX_56" # Participant id of interest
markers <- names(fs_VHX_BPD) # Names of markers
n_markers <- length(markers) # Number of markers
# All possible permutations:
all_permutations <- gtools::permutations(n = n_markers, r = n_markers)
n_all <- nrow(all_permutations) # Number of all possible permutations
writeLines(sprintf("There are %s possible marker permutations", n_all))
n_plot <- 5 # Number of permutations to plot 
inds <- sample(x = n_all, size = n_plot, replace = F) # Indices of permutations

results <- lapply(inds, function(ind){ # For n_plot permutations of markers
  permuted_markers <- markers[all_permutations[ind,]] # for permuted markers
  results_m <- lapply(1:length(markers), function(m){ # for marker subsets
    # extract data on maker subset
    y <- lapply(ys_VHX_BPD[[pid]], function(y_epi) y_epi[permuted_markers][1:m]) 
    # compute posterior probabilities
    z <- compute_posterior(y, fs = fs_VHX_BPD, return.RG = TRUE, return.logp = TRUE)
    # extract log likelihoods of relationship graphs RG
    llikes <- sapply(z$RGs, function(RG) RG$logp)
    # extract maximum likelihood (ML) RG 
    RG <- z$RGs[[which(abs(llikes - max(llikes)) < .Machine$double.eps^0.5)[1]]]
    # return probabilities and ML RG
    list(probs = z$marg, RG = RG_to_igraph(RG, determine_MOIs(y)))
  })
  # extract and return probs and RGs
  relapse_probs <- sapply(results_m, function(x) x[["probs"]])[2,] 
  MLRGs <- sapply(results_m, function(x) x[["RG"]]) 
  list(relapse_probs = relapse_probs, MLRGs = MLRGs)
})

#===============================================================================
# Plot data
#===============================================================================
oldpar <- par(no.readonly = T) # To restore plotting parameters later
plot_data(ys_VHX_BPD[pid], fs = fs_VHX_BPD) # Plot Data

#===============================================================================
# Plot evolution of reinfection probability over n_markers
#===============================================================================
cols <- RColorBrewer::brewer.pal(min(8, n_plot), "Dark2")
par(mar = c(4,5,1,2), pty = "m")
plot(NULL, ylim = c(0,1), xlim = c(1,n_markers), bty = "n",
     xlab = "Number of markers", ylab = "Probability of relapse")
for(i in 1:n_plot){
  lines(y = results[[i]][["relapse_probs"]], x = 1:n_markers, 
       col = cols[i], type = "b", pch = 20)
}

################################################################################
# Plot ML RGs (not visually useful when n_markers is large)
################################################################################
par(mar = rep(.5,4), pty = "s", mfrow = c(n_plot,n_markers))
for(i in 1:n_plot){
  for(j in 1:n_markers){
  plot_RG(results[[i]][["MLRGs"]][[j]], vertex.size = 20, 
          vertex.label = NA, vertex.palette = "Greys")
  title(main = markers[all_permutations[inds[i],]][j], cex.main = 0.5) 
  box(col = cols[i])
  }
}

par(oldpar)
