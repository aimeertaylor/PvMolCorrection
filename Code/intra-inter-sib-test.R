################################################################################
# Write a function to test if cases of suspecting misspecification align with 
# cases where RGS with sibling edges both within and across infections have -INF
# likelihood
################################################################################
library(Pv3Rs)
rm(list = ls())
y <- ys_VHX_BPD[["VHX_56"]]
plot_data(ys_VHX_BPD["VHX_56"], fs = fs_VHX_BPD) # Data suggest possible relapse
post <- compute_posterior(y, fs = fs_VHX_BPD, return.RG = TRUE, return.logp = TRUE)
post$marg # Low probability of reinfection

# Explore likelihoods: 
# Extract all log likelihoods in decreasing order
llikes <- (sapply(post$RGs, function(RG) RG$logp))
sorted_llikes <- sort(llikes, decreasing = T) # Sort log likelihoods
plot(sorted_llikes[!is.infinite(sorted_llikes)]) # Plot sorted finite log likelihoods

# Plot RGs with maximum log likelihoods
for(i in 1:length(unique(sorted_llikes[!is.infinite(sorted_llikes)]))){
  RGs <- post$RGs[which(abs(llikes - unique(sorted_llikes)[i]) < .Machine$double.eps^0.5)]
  par(mar = rep(0.1,4), mfrow = c(1,length(RGs)))
  for(j in 1:length(RGs)) {
    plot_RG(RG_to_igraph(RGs[[j]], determine_MOIs(y)), vertex.size = 20)
    box()
  }  
}

# Write a function that explores the absence of within both and across


# Plot RGs with next largest log likelihood
RGs <- post$RGs[which(abs(llikes - unique(sorted_llikes)[2]) < .Machine$double.eps^0.5)]
par(mar = rep(0.1,4), mfrow = c(1,length(RGs)))
for(i in 1:length(RGs)) {
  plot_RG(RG_to_igraph(RGs[[i]], determine_MOIs(y)), vertex.size = 20)
  box()
}

# Plot RGs with next largest log likelihood
RGs <- post$RGs[which(abs(llikes - unique(sorted_llikes)[3]) < .Machine$double.eps^0.5)]
par(mar = rep(0.1,4), mfrow = c(1,length(RGs)))
for(i in 1:length(RGs)) {
  plot_RG(RG_to_igraph(RGs[[i]], determine_MOIs(y)), vertex.size = 20)
  box()
}


# Now let's explore the equivalence class with the largest log likelihood.
# In the following code, we place two graphs in the same equivalence class if
# they share the same likelihood. This is not ideal (two graphs that are not
# isomorphic up to permutation could share the same likelihood), but it works
# here: the plot shows only isomorphic graphs within the equivalence class.
adj_equal <- abs(diff(sorted_llikes, lag = 1)) < .Machine$double.eps^0.5 # Find matches
decr_idxs <- which(adj_equal == FALSE) # Change points: 3, 9, 15, ...
class_sizes <- c(decr_idxs[1], diff(decr_idxs)) # Number of graphs per class

# log likelihood of representative from each 'equivalence class' (EC)
llikes_unique <- sorted_llikes[decr_idxs]

# EC likelihood
class_ps <- exp(llikes_unique)*class_sizes
max_class_p <- which(class_ps == max(class_ps)) # ML EC index 
max_idx <- decr_idxs[max_class_p] # Index of last graph in ML EC
max_size <- class_sizes[max_class_p] # Number of graphs in ML EC

# Plot all graphs within the ML EC 
par(mar = rep(0.1,4), mfrow = c(1,3))
RG_order <- order(llikes, decreasing = T) # order RGs by logl
for(i in (max_idx-max_size+1):max_idx) { # EC consists of the RGs with logl rank 21-32
  RG <- post$RGs[[RG_order[i]]]
  RG_igraph <- RG_to_igraph(RG, determine_MOIs(y))
  plot_RG(RG_igraph, vertex.size = 25, vertex.label = NA) 
  box()
}

