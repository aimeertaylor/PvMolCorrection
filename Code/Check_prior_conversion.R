####### This is a mess

# check prior conversion, section 2.2, supplement of Foo et al Bioinformatics
rm(list = ls())
load("../RData/prior_estimates.RData")
load("../RData/marg_results_Pv3Rs.RData")
load("../RData/results_Pv3Rs.RData")
pids <- names(ys_VHX_BPD)
n_episodes <- sapply(ys_VHX_BPD, length)
states <- c("C", "L", "I")

# Check episodes are ordered correctly
all(sapply(ys_VHX_BPD, function(x) all(sort(as.numeric(names(x))) == as.numeric(names(x)))))

approx <- t(sapply(pids, function(pid){
  #*** 
  pid <- pids[100]
  
  eids <- paste(pid, names(ys_VHX_BPD[[pid]]), sep = "_")
  rids <- eids[-1] # recurrences ids (works only if check above passes)
  sequences <- gtools::permutations(n = length(states), v = states, r = n_episodes[pid]-1, repeats.allowed = T)
  marg_prior <- prior[rids,]
  marg_post <- Uniform_Pv3Rs[rids,]
  
  approx_unnorm <- apply(sequences, 1, function(sequence){
    seq_prior <- marg_prior[1, sequence[1]] * marg_prior[2, sequence[2]] * marg_prior[3, sequence[3]]
    seq_post <- marg_post[1, sequence[1]] * marg_post[2, sequence[2]] * marg_post[3, sequence[3]]
    seq_post
  })
  
  sequence_probs <- approx_unnorm / sum(approx_unnorm)
  names(sequence_probs) <- apply(sequences, 1, paste, collapse = '')
  plot(Uniform_joint[[pid]]$joint[names(sequence_probs)], sequence_probs)
  abline(a = 0, b = 1)
  
  c("C" = sum(sequences[sequences[,1] == "C",4]), 
    "L" = sum(sequences[sequences[,1] == "L",4]), 
    "I" = sum(sequences[sequences[,1] == "I",4]))
    
  TimeToEvent_Pv3Rs[rids,]
    
    
  approx_unnorm <- Uniform_Pv3Rs[reids, states] * prior[eid, states]
  approx_unnorm / sum(approx_unnorm)
}))

# Almost perfect
par(mfrow = c(2,2))
for(s in states){
  plot(y = approx[eids,s], x = TimeToEvent_Pv3Rs[eids, s], pch = 20, cex = 0.5)
}

