rm(list = ls())
source("match_counting.R")
source("genetic_proximity.R")

# Compute number of episodes per participant 
n_epi_per_pid <- sapply(ys_VHX_BPD, length)

# Compute matches for pids with one or more recurrences 
matches <- sapply(ys_VHX_BPD[n_epi_per_pid > 1], match_counting)

# Compute genetic proximity for pids with one or more recurrences
proxies_one <- sapply(ys_VHX_BPD[n_epi_per_pid > 1], genetic_proximity, type = "rhat_one")
proxies_avg <- sapply(ys_VHX_BPD[n_epi_per_pid > 1], genetic_proximity, type = "rhat_avg")
proxies_tot <- sapply(ys_VHX_BPD[n_epi_per_pid > 1], genetic_proximity, type = "rhat_tot")

# Extract values into vector form
matches_p <- unlist(sapply(matches, function(x) x["p",])) # match proportion
rhats_one <- unlist(sapply(proxies_one, function(x) x["rhat",])) # rhat
rhats_avg <- unlist(sapply(proxies_avg, function(x) x["rhat",])) # rhat
rhats_tot <- unlist(sapply(proxies_tot, function(x) x["rhat",])) # rhat

matches_n <- unlist(sapply(matches, function(x) x["n",])) # match marker count
rhats_one_n <- unlist(sapply(proxies_one, function(x) x["nmarks",])) # rhat marker count
rhats_avg_n <- unlist(sapply(proxies_avg, function(x) x["nmarks",])) # rhat marker count
rhats_tot_n <- unlist(sapply(proxies_tot, function(x) x["nmarks",])) # rhat marker count

# Name extracted values; n.b., gsub(".", "_") does not work because pids with
# one recurrence miss "_2" e.g. BPD_103
names(matches_p) <- names(matches_n) <- unlist(sapply(names(matches), function(pid) {
  paste(pid, colnames(matches[[pid]]), sep = "_") 
}))
names(rhats_one) <- names(rhats_one_n) <- unlist(sapply(names(proxies_one), function(pid) {
  paste(pid, colnames(proxies_one[[pid]]), sep = "_") 
}))
names(rhats_avg) <- names(rhats_avg_n) <- unlist(sapply(names(proxies_avg), function(pid) {
  paste(pid, colnames(proxies_avg[[pid]]), sep = "_") 
}))
names(rhats_tot) <- names(rhats_tot_n) <- unlist(sapply(names(proxies_tot), function(pid) {
  paste(pid, colnames(proxies_tot[[pid]]), sep = "_") 
}))

# Save extracted values
save(matches_p, rhats_one, rhats_avg, rhats_tot, 
     matches_n, rhats_one_n, rhats_avg_n, rhats_tot_n, file = "../RData/genetic_proximities.RData")
