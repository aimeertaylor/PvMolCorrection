library(Pv3Rs)
source("match_counting.R")

# Compute number of episodes per participant 
n_epi_per_pid <- sapply(ys_VHX_BPD, function(x) length(x))

# Compute matches for  paired data 
matches_n_p <- sapply(ys_VHX_BPD[n_epi_per_pid > 1], match_counting)
matches_p <- unlist(sapply(matches_n_p, function(x) x["p",]))
matches_n <- unlist(sapply(matches_n_p, function(x) x["n",]))

# I first did names(matches_n) <- gsub(".", "_", names(matches_n), fixed = T)
# Need following because otherwise pids with one episode miss "_2" e.g. BPD_103
names(matches_p) <- names(matches_n) <- unlist(sapply(names(matches_n_p), function(n) {
  paste(n, colnames(matches_n_p[[n]]), sep = "_")
}))

# Load posterior probabilities
load("../RData/marg_results_Pv3Rs.RData")

if (!all(names(matches_n) == row.names(Uniform_Pv3Rs))) stop()

# Plot match proportion against 
plot(x = matches_p, y = Uniform_Pv3Rs[,"L"], 
     cex = matches_n/9,
     xlab = "Match proportion", 
     ylab = "Posterior probability (based on genetic data alone)")

plot(x = matches_p, y = Uniform_Pv3Rs[,"C"], 
     cex = matches_n/9,
     xlab = "Match proportion", 
     ylab = "Posterior probability (based on genetic data alone)")

plot(x = matches_p, y = Uniform_Pv3Rs[,"I"], # Use this plot for diagnosis
     cex = matches_n/9,
     xlab = "Match proportion", 
     ylab = "Posterior probability (based on genetic data alone)")
# Annotate already identified cases

half_sibs <- c("VHX_39_2", "VHX_56_2", "VHX_91_2", "VHX_113_6", "VHX_329_4", "VHX_450_8", 
               "VHX_489_4", "VHX_529_4", "VHX_532_4", "VHX_225_4")

text(x = matches_p[half_sibs], y = Uniform_Pv3Rs[half_sibs,"I"], labels = half_sibs)
abline(h = 0.5, v = 4/9, lty = "dotted")
abline(v = (0:9)/9, lty = "dotted", col = "red")

suspect_epis <- names(which(matches_p >= 4/9 & Uniform_Pv3Rs[,"I"] > 0.5))
suspect_pids <- unique(apply(do.call(rbind, strsplit(suspect_epis, split = "_"))[,1:2], 1, paste, collapse = "_"))
for(suspect_pid in suspect_pids) {
  plot_data(ys = ys_VHX_BPD[suspect_pid], fs = fs_VHX_BPD)
  mtext(paste(suspect_epis[grepl(suspect_pid, suspect_epis)], 
        collapse = " "))
}

