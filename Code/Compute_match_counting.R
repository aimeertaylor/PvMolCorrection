library(Pv3Rs)
source("match_counting.R")
source("genetic_proximity.R")
load("../RData/marg_results_Pv3Rs.RData") # Load posterior probabilities

# Compute number of episodes per participant 
n_epi_per_pid <- sapply(ys_VHX_BPD, length)

# Compute matches for pids with one or more recurrences 
matches <- sapply(ys_VHX_BPD[n_epi_per_pid > 1], match_counting)

# Compute genetic proximity for pids with one or more recurrences
proxies_one <- sapply(ys_VHX_BPD[n_epi_per_pid > 1], genetic_proximity, type = "rhat_one")
proxies_avg <- sapply(ys_VHX_BPD[n_epi_per_pid > 1], genetic_proximity, type = "rhat_avg")

# Extract values into vector form
matches_p <- unlist(sapply(matches, function(x) x["p",])) # match proportion
rhats_one <- unlist(sapply(proxies_one, function(x) x["rhat",])) # rhat
rhats_avg <- unlist(sapply(proxies_avg, function(x) x["rhat",])) # rhat

matches_n <- unlist(sapply(matches, function(x) x["n",])) # match marker count
rhats_one_n <- unlist(sapply(proxies_one, function(x) x["nmarks",])) # rhat marker count
rhats_avg_n <- unlist(sapply(proxies_avg, function(x) x["nmarks",])) # rhat marker count


# Name extracted values; n.b., gsub(".", "_") does not work because pids with
# one recurrence miss "_2" e.g. BPD_103
names(matches_p) <- names(matches_n) <- unlist(sapply(names(matches), function(pid) {
  paste(pid, colnames(matches_n_p[[pid]]), sep = "_") 
}))
names(rhats_one) <- names(rhats_one_n) <- unlist(sapply(names(proxies_one), function(pid) {
  paste(pid, colnames(proxies_one[[pid]]), sep = "_") 
}))
names(rhats_avg) <- names(rhats_avg_n) <- unlist(sapply(names(proxies_avg), function(pid) {
  paste(pid, colnames(proxies_avg[[pid]]), sep = "_") 
}))


if (!all(names(matches_p) == row.names(Uniform_Pv3Rs))) stop()
if (!all(names(rhats_one) == row.names(Uniform_Pv3Rs))) stop()
if (!all(names(rhats_avg) == row.names(Uniform_Pv3Rs))) stop()

half_sibs <- c("VHX_39_2", "VHX_56_2", "VHX_91_2", "VHX_113_6", "VHX_450_8", 
               "VHX_489_4", "VHX_529_4", "VHX_532_4", "VHX_225_4")

pdf("Pv3Rs vs distance.pdf")
par(mfrow = c(2,2), pty = "s", height = 12, width = 12)

# # Plot match proportion 
plot(x = matches_p, y = 1-Uniform_Pv3Rs[,"I"],
     cex = matches_n/9, pch = 19, bty = "n",
     xlab = "genetic proximity: match proportion",
     ylab = "Posterior relapse plus recrudescence probability (uniform prior on states)")
text(x = matches_p[half_sibs], y = 1-Uniform_Pv3Rs[half_sibs,"I"],
     labels = half_sibs, cex = 0.5, pos = 3)
abline(h = 0.5, v = 4/9, lty = "dashed")
abline(v = (0:9)/9, lty = "dotted", col = "grey")

# Plot genetic proximity assuming one related genotype
plot(x = rhats_one, y = 1-Uniform_Pv3Rs[,"I"], 
     cex = rhats_one_n/9, pch = 19, bty = "n", 
     xlab = "genetic proximity: relatedness assuming one related genotype", 
     ylab = "Posterior relapse plus recrudescence probability (uniform prior on states)")
text(x = rhats_one[half_sibs], y = 1-Uniform_Pv3Rs[half_sibs,"I"], 
     labels = half_sibs, cex = 0.5, pos = 3)
abline(h = 0.5, v = 0.5, lty = "dashed")
# Add additional suspect 
suspect_epi_one <- names(which(rhats_one >= 0.5 & Uniform_Pv3Rs[,"I"] > 0.5))
extra <- suspect_epi_one[!suspect_epi_one %in% half_sibs]
text(x = rhats_one[extra], y = 1-Uniform_Pv3Rs[extra,"I"], 
     labels = extra, cex = 0.5, pos = 3, col = "red")


# Plot genetic proximity: average relatedness
plot(x = rhats_avg, y = 1-Uniform_Pv3Rs[,"I"], 
     cex = rhats_avg_n/9, pch = 19, bty = "n", 
     xlab = "genetic proximity: relatedness averaged over genotypes", 
     ylab = "Posterior relapse plus recrudescence probability (uniform prior on states)")
text(x = rhats_avg[half_sibs], y = 1-Uniform_Pv3Rs[half_sibs,"I"], 
     labels = half_sibs, cex = 0.5, pos = 3)
abline(h = 0.5, v = 0.5, lty = "dashed")
# Add additional suspect 
suspect_epi_avg <- names(which(rhats_avg >= 0.5 & Uniform_Pv3Rs[,"I"] > 0.5))
extra <- suspect_epi_avg[!suspect_epi_avg %in% half_sibs]
text(x = rhats_avg[extra], y = 1-Uniform_Pv3Rs[extra,"I"], 
     labels = extra, cex = 0.5, pos = 3, col = "red")
dev.off()


# Inspect data 
all(suspect_epi_avg %in% suspect_epi_one)
suspect_pids_one <- unique(apply(do.call(rbind, strsplit(suspect_epi_one, split = "_"))[,1:2], 1, paste, collapse = "_"))
plot_data(ys = ys_VHX_BPD[suspect_pids_one], fs = fs_VHX_BPD)

