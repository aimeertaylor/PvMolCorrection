library(Pv3Rs)
source("match_counting.R")

# Compute number of episodes per participant 
n_epi_per_pid <- sapply(ys_VHX_BPD, function(x) length(x))

# Compute matches for  paired data 
matches_n_p <- sapply(ys_VHX_BPD[n_epi_per_pid > 1], match_counting)
matches_n <- unlist(sapply(matches_n_p, function(x) x["p",]))

# I first did names(matches_n) <- gsub(".", "_", names(matches_n), fixed = T)
# Need following because otherwise pids with one episode miss "_2" e.g. BPD_103
names(matches_n) <- unlist(sapply(names(matches_n_p), function(n) {
  paste(n, colnames(matches_n_p[[n]]), sep = "_")
}))

# Load posterior probabilities
load("../RData/marg_results_Pv3Rs.RData")

if (!all(names(matches_n) == row.names(Uniform_Pv3Rs))) stop()


# Plot match proportion against 
plot(x = matches_n, y = Uniform_Pv3Rs[,"L"], 
     xlab = "Match proportion", 
     ylab = "Posterior relapse probability (based on genetic data alone)")

plot(x = matches_n, y = Uniform_Pv3Rs[,"C"], 
     xlab = "Match proportion", 
     ylab = "Posterior relapse probability (based on genetic data alone)")

plot(x = matches_n, y = Uniform_Pv3Rs[,"I"], 
     xlab = "Match proportion", 
     ylab = "Posterior relapse probability (based on genetic data alone)")

# Something wrong? 
suspect_epis <- names(which(matches_n == 0 & Uniform_Pv3Rs[,"L"] > 0.5))
suspect_pids <- unique(apply(do.call(rbind, strsplit(suspect_epis, split = "_"))[,1:2], 1, paste, collapse = "_"))

plot_data(ys = ys_VHX_BPD[suspect_pids], fs = fs_VHX_BPD)
