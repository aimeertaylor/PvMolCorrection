################################################################################
# Prior conversion but with access to only marginal probabilities
################################################################################
rm(list = ls())
load("../RData/prior_estimates.RData")
load("../RData/marg_results_Pv3Rs.RData")
states <- c("C", "L", "I")
n_episodes <- sapply(ys_VHX_BPD, length)
pids_plus <- names(which(n_episodes > 2))
eids <- rownames(Uniform_Pv3Rs)

# Prior conversion using equation in section 2.2, supplement of Foo et al 2026
TimeToEvent_approx <- t(sapply(eids, function(eid){
  approx_unnorm <- Uniform_Pv3Rs[eid, states] * prior[eid, states]
  approx_unnorm / sum(approx_unnorm)
}))

# Ascertain which recurrences are in individuals where marginal distribution is not joint
eids_plus <- sapply(eids, function(eid) {
  pid_eid <- paste(strsplit(eid, split = "_")[[1]][1:2], collapse ="_")
  pid_eid %in% pids_plus
})

# Almost perfect prior conversion when marginal == joint
par(mfrow = c(2,2))
for(s in states){
  plot(y = TimeToEvent_approx[!eids_plus,s], 
       x = TimeToEvent_Pv3Rs[!eids_plus, s], 
       ylab = "Marginal probability approximation", 
       xlab = "Jointly modelled probability",
       xlim = c(0,1), ylim = c(0,1),
       pch = 20, cex = 0.5, bty = "n", main = s)
}

# Difference is tiny (same order as machine precision)
range(unlist(unlist(TimeToEvent_approx[!eids_plus,states])) - 
        unlist(as.vector(TimeToEvent_Pv3Rs[!eids_plus,states])))

# Imperfect prior conversion when marginal != joint
par(mfrow = c(2,2))
for(s in states){
  plot(y = TimeToEvent_approx[eids_plus,s], 
       x = TimeToEvent_Pv3Rs[eids_plus, s], 
       ylab = "Marginal probability approximation", 
       xlab = "Jointly modelled probability",
       xlim = c(0,1), ylim = c(0,1),
       pch = 20, cex = 0.5, bty = "n", main = s)
}

# Difference is small: 
range(unlist(unlist(TimeToEvent_approx[eids_plus,states])) - 
        unlist(as.vector(TimeToEvent_Pv3Rs[eids_plus,states])))
