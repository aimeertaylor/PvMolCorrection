rm(list = ls())
library(Pv3Rs)
source("match_counting.R")
source("genetic_proximity.R")
load("../RData/quantiles_null.RData")
load("../RData/marg_results_Pv3Rs.RData") # Load posterior probabilities
Figs <- T

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

# Check
if (!all(names(matches_p) == row.names(Uniform_Pv3Rs))) stop()
if (!all(names(rhats_one) == row.names(Uniform_Pv3Rs))) stop()
if (!all(names(rhats_avg) == row.names(Uniform_Pv3Rs))) stop()
if (!all(names(rhats_tot) == row.names(Uniform_Pv3Rs))) stop()

# Plot probabilities against genetic proximity: 
par(mfrow = c(2,2), pty = "s")

# Match proportion 
plot(x = matches_p, y = 1-Uniform_Pv3Rs[,"I"],
     cex = matches_n/9, pch = 19, bty = "n",
     xlab = "Match proportion",
     ylab = "Relapse plus recrudescence probability")
abline(h = q_3Rs, v = q_mat, lty = "dashed", col = "red")
abline(v = (0:9)/9, lty = "dotted", col = "grey")

# Average relatedness
plot(x = rhats_avg, y = 1-Uniform_Pv3Rs[,"I"], 
     cex = rhats_avg_n/9, pch = 19, bty = "n", 
     xlab = "Relatedness averaged over genotypes", 
     ylab = "Relapse plus recrudescence probability")
abline(h = q_3Rs, v = q_avg, lty = "dashed", col = "red")

# Assuming one related genotype
plot(x = rhats_one, y = 1-Uniform_Pv3Rs[,"I"], 
     cex = rhats_one_n/9, pch = 19, bty = "n", 
     xlab = "Relatedness assuming one related genotype", 
     ylab = "Relapse plus recrudescence probability")
abline(h = q_3Rs, v = q_one, lty = "dashed", col = "red")

# Total relatedness
plot(x = rhats_tot, y = 1-Uniform_Pv3Rs[,"I"], 
     cex = rhats_one_n/9, pch = 19, bty = "n", 
     xlab = "Total relatedness", 
     ylab = "Relapse plus recrudescence probability")
abline(h = q_3Rs, v = q_tot, lty = "dashed", col = "red")

#===============================================================================
# Which distance identifies most potential FN: tot and one and one is prettier
#===============================================================================
FN_tot <- names(which(rhats_tot > q_tot & (1-Uniform_Pv3Rs[,"I"]) < q_3Rs))
FN_one <- names(which(rhats_one > q_one & (1-Uniform_Pv3Rs[,"I"]) < q_3Rs))
FN_avg <- names(which(rhats_avg > q_avg & (1-Uniform_Pv3Rs[,"I"]) < q_3Rs))

length(FN_tot) # 24
length(FN_one) # 24
length(FN_avg) # 23

all(FN_tot %in% FN_one) & all(FN_avg %in% FN_one)

#===============================================================================
# Plot for ms
#===============================================================================
par(mfrow = c(1,2))
pids <- apply(do.call(rbind, strsplit(names(rhats_one), split = "_"))[,1:2], 1, paste, collapse = "_")
# Group genetic data by patient ID:
pid_arm <- unique(do.call(rbind, plyr::dlply(Combined_Time_Data[, c("patientid", "arm_num")], "patientid")))
PMQ <- grepl("PMQ", pid_arm[,"arm_num"])
names(PMQ) <- pid_arm[,"patientid"]

plot(x = rhats_one, y = 1-Uniform_Pv3Rs[,"I"], 
     cex = rhats_one_n/9, pch = 19, bty = "n", 
     col = PMQ[pids] + 1, 
     xlab = "Genetic proximity", 
     ylab = "Relapse plus recrudescence probability (uniform prior)")
abline(v = q_one, lty = "dashed", col = "red")
text(x = rhats_one[FN_one], y = 1-Uniform_Pv3Rs[FN_one,"I"], 
     labels = FN_one, pos = 3, cex = 0.5)

plot(x = rhats_one, y = 1-TimeToEvent_Pv3Rs[,"I"], 
     cex = rhats_one_n/9, pch = 19, bty = "n", 
     col = PMQ[pids] + 1, 
     xlab = "Genetic proximity: relatedness assuming one related genotype", 
     ylab = "Relapse plus recrudescence probability (informative prior)")
abline(v = q_one, lty = "dashed", col = "red")
text(x = rhats_one[FN_one], y = 1-TimeToEvent_Pv3Rs[FN_one,"I"], 
     labels = FN_one, pos = 3, cex = 0.5)

#===============================================================================
# Numbers of ms: 
#===============================================================================
# Average probability of relapse among reinfection-classified 
sum(Uniform_Pv3Rs[which(rhats_one <= q_one), "L"])/sum(rhats_one <= q_one)
sum(TimeToEvent_Pv3Rs[which(rhats_one <= q_one), "L"])/sum(rhats_one <= q_one)

# Aside: average probability of recrudescence among reinfection-classified 
sum(Uniform_Pv3Rs[which(rhats_one <= q_one), "C"])/sum(rhats_one <= q_one)
sum(TimeToEvent_Pv3Rs[which(rhats_one <= q_one), "C"])/sum(rhats_one <= q_one)

# Proportion / probability of "relapse or recrudescence" recurrences among
# recurrences in PMQ treated participants:
sum((1-Uniform_Pv3Rs[,"I"])*PMQ[pids])/sum(PMQ[pids])
sum((1-TimeToEvent_Pv3Rs[,"I"])*PMQ[pids])/sum(PMQ[pids])
sum(rhats_one*PMQ[pids])/sum(PMQ[pids])

# Proportion / probability of "relapse or recrudescence" recurrences among
# reinfection-classified recurrences in PMQ treated participants:
sum((1-Uniform_Pv3Rs[,"I"])*PMQ[pids]*(rhats_one <= q_one))/sum(PMQ[pids]*(rhats_one <= q_one))
sum((1-TimeToEvent_Pv3Rs[,"I"])*PMQ[pids]*(rhats_one <= q_one))/sum(PMQ[pids])

#===============================================================================
# Inspect data on outliers: 
#===============================================================================
pids <- unique(apply(do.call(rbind, strsplit(FN_one, split = "_"))[,1:2], 1, paste, collapse = "_")) 

# Extract all episodes for pids 
episodes <- unname(unlist(sapply(pids, function(pid) {
  sapply(names(ys_VHX_BPD[[pid]]), function(epi) {
    paste(pid, epi, sep = "_")})})))

# Extract probabilities for all recurrences in pids
all_recs_pids <- unlist(sapply(pids, function(pid) paste(pid, names(ys_VHX_BPD[[pid]])[-1], sep = "_")))
LorC <- 1-Uniform_Pv3Rs[all_recs_pids,"I"]
names(LorC) <- all_recs_pids 
no_episodes_per_pid <- sapply(pids, function(pid) length(ys_VHX_BPD[[pid]]))
text <- rep("", length(episodes)) 
names(text) <- episodes
text[names(LorC)] <- round(LorC, 2)*100
text_col <- unlist(sapply(1:length(pids), function(i){
  if(i %% 2 == 1) {
    rep("white", no_episodes_per_pid[i])
  } else {
    rep("black", no_episodes_per_pid[i])
  }}))
main_mar <- c(3, 3.5, 1.5, 3.5)
z <- 1/(2*length(episodes)) # See text x placement

# Plot for ms
if(Figs) png("../Figures/data_false_neg_failures.png", width = 10, height = 7, units = "in", res = 300) 
par(fig = c(0,1,0.2+0.01,1), mar = main_mar, mfrow = c(1,1)) # Important to call before and after plot_data
plot_data(ys = ys_VHX_BPD[pids], fs = fs_VHX_BPD, marker.annotate = F, mar = main_mar, person.vert = T)
par(fig = c(0,1,0.2+0.01,1), mar = main_mar) # Important to call before and after plot_data
text(y = rep(0.0375, length(episodes)),
     x = seq(z, 1-z, length.out = length(episodes)), 
     labels = text, cex = 0.4, col = text_col, srt = 90)
points(y = rep(0.06, length(unlist(FN_one))),
       x = seq(z, 1-z, length.out = length(episodes))[episodes %in% unlist(FN_one)], 
       pch = 25, cex = 0.4, bg = "black") 
if(Figs) dev.off()


