rm(list = ls())
library(Pv3Rs)
load("../RData/genetic_proximities.RData")
load("../RData/quantiles_null.RData")
load("../RData/marg_results_Pv3Rs.RData") # Load posterior probabilities
Figs <- T

#===============================================================================
# Plot probabilities against genetic proximity to decide on which proximity: 
#===============================================================================
par(mfrow = c(2,2), pty = "s")

# Match proportion 
plot(x = matches_p, y = 1-Uniform_Pv3Rs[,"I"],
     cex = matches_n/9, pch = 19, bty = "n",
     xlab = "Match proportion",
     ylab = "Relapse plus recrudescence probability")
abline(h = q_3Rs, v = q_mat, lty = c("dashed", "dotted"), col = "red")
abline(v = (0:9)/9, lty = "dotted", col = "grey")

# Average relatedness
plot(x = rhats_avg, y = 1-Uniform_Pv3Rs[,"I"], 
     cex = rhats_avg_n/9, pch = 19, bty = "n", 
     xlab = "Relatedness averaged over genotypes", 
     ylab = "Relapse plus recrudescence probability")
abline(h = q_3Rs, v = q_avg, lty = c("dashed", "dotted"), col = "red")

# Assuming one related genotype
plot(x = rhats_one, y = 1-Uniform_Pv3Rs[,"I"], 
     cex = rhats_one_n/9, pch = 19, bty = "n", 
     xlab = "Relatedness assuming one related genotype", 
     ylab = "Relapse plus recrudescence probability")
abline(h = q_3Rs, v = q_one, lty = c("dashed", "dotted"), col = "red")

# Total relatedness
plot(x = rhats_tot, y = 1-Uniform_Pv3Rs[,"I"], 
     cex = rhats_tot_n/9, pch = 19, bty = "n", 
     xlab = "Total relatedness", 
     ylab = "Relapse plus recrudescence probability")
abline(h = q_3Rs, v = q_tot, lty = c("dashed", "dotted"), col = "red")

#===============================================================================
# Which distance identifies most potential FN at 95% significance: tot and one
# Use total relatedness hereafter
#===============================================================================
FN_tot <- names(which(rhats_tot > q_tot["95%"] & (1-Uniform_Pv3Rs[,"I"]) < q_3Rs["95%"]))
FN_one <- names(which(rhats_one > q_one["95%"] & (1-Uniform_Pv3Rs[,"I"]) < q_3Rs["95%"]))
FN_avg <- names(which(rhats_avg > q_avg["95%"] & (1-Uniform_Pv3Rs[,"I"]) < q_3Rs["95%"]))

length(FN_tot) # 24
length(FN_one) # 24
length(FN_avg) # 23

all(FN_one %in% FN_tot) & all(FN_avg %in% FN_tot) # total relatedness covers all


#===============================================================================
# Load data from all participants to extract treatment
#===============================================================================
load("../jwatowatson-RecurrentVivax-4870715/RData/TimingModel/Combined_Time_Event.RData") 
PMQ <- grepl("PMQ", Combined_Time_Data$arm_num)
names(PMQ) <- paste(Combined_Time_Data$patientid, Combined_Time_Data$episode, sep = "_")

#===============================================================================
# Load and inspect names of discrepant recurrence probabilities
#===============================================================================
load("../RData/big_diffs_pv3rs_prototype.RData")
load("../RData/big_diffs_joint_pairwise.RData")
eids_proto_pv3Rs
eids_joint_pwise

#===============================================================================
# Load and extract Pv3Rs pairwise
#===============================================================================
load("../RData/results_Pv3Rs.RData")
Uniform_pwise <- do.call(rbind, Uniform_pairwise)
rownames(Uniform_pwise) <- unlist(sapply(names(Uniform_pairwise), function(pid){
  paste(pid, rownames(Uniform_pairwise[[pid]]), sep = "_")
}))


#===============================================================================
# Compare Pv3Rs (data modelled jointly and pairwise) with prototype
#
# Remember that probabilities given uniform prior were only computed using the
# prototype for a subset of recurrences. The subset does not include
# eids_joint_pwise - see below. For all recurrences in the subset, Pv3Rs was
# able to model data jointly; also see below.
#
# Take away: 
# Pv3Rs generally smoother 
# Pv3Rs improves VHX_452_2 and VHX_214_2
# Pv3Rs clearly worse for VHX_56_2, VHX_91_2 
# Pv3Rs possibly worse for BPD_70_2, VHX_298_2, BPD_253_3  
# Pv3Rs pairwise: minor improvement for BPD_253_3  
#===============================================================================
load('../jwatowatson-RecurrentVivax-4870715/RData/GeneticModel/Including_Complex_Cases_Full_Posterior_Model_samples_Tagnostic.RData')
any(eids_joint_pwise %in% rownames(thetas_9MS_Tagnostic)) # subset does not inc. eids_joint_pwise
uls <- grepl("Uniform", names(eids_proto_pv3Rs)) # Logical of discrepant uniform
par(mfrow = c(3,1), pty = "m")

# Prototype
plot(x = rhats_tot[rownames(thetas_9MS_Tagnostic)], 
     y = 1-thetas_9MS_Tagnostic[,"I"], 
     pch = 19, bty = "n", 
     cex = rhats_tot_n[rownames(thetas_9MS_Tagnostic)]/9,
     xlab = "Genetic proximity: total relatedness", 
     ylab = "Prototype")
abline(v = q_tot, lty = c("dashed", "dotted"), col = "red")
text(x = rhats_tot[eids_proto_pv3Rs[uls]], 
     y = 1-thetas_9MS_Tagnostic[eids_proto_pv3Rs[uls],"I"], 
     labels = eids_proto_pv3Rs[uls], pos = 4, cex = 0.5, col = "hotpink")

# Pv3Rs comparable recurrences computed jointly:
all(Uniform_Pv3Rs[rownames(thetas_9MS_Tagnostic), "joint"]) 
plot(x = rhats_tot[rownames(thetas_9MS_Tagnostic)], 
     y = 1-Uniform_Pv3Rs[rownames(thetas_9MS_Tagnostic),"I"], 
     pch = 19, bty = "n", 
     cex = rhats_tot_n[rownames(thetas_9MS_Tagnostic)]/9,
     xlab = "Genetic proximity: total relatedness", 
     ylab = "Pv3Rs (data modelled jointly)")
abline(v = q_tot, lty = c("dashed", "dotted"), col = "red")
text(x = rhats_tot[eids_proto_pv3Rs[uls]], 
     y = 1-Uniform_Pv3Rs[eids_proto_pv3Rs[uls],"I"], 
     labels = eids_proto_pv3Rs[uls], pos = 4, cex = 0.5, col = "hotpink")

# Pv3Rs comparable recurrences computed pairwise
plot(x = rhats_tot[rownames(thetas_9MS_Tagnostic)], 
     y = 1-Uniform_pwise[rownames(thetas_9MS_Tagnostic), "I"],
     pch = 19, bty = "n", 
     cex = rhats_tot_n[rownames(thetas_9MS_Tagnostic)]/9, 
     xlab = "Genetic proximity: total relatedness", 
     ylab = "Pv3Rs (data modelled pairwise)")
abline(v = q_tot, lty = c("dashed", "dotted"), col = "red")
text(x = rhats_tot[eids_proto_pv3Rs[uls]], y = 1-Uniform_pwise[eids_proto_pv3Rs[uls],"I"], 
     labels = eids_proto_pv3Rs[uls], pos = 4, cex = 0.5, col = "hotpink")

#===============================================================================
# Plots for ms
#===============================================================================
# List recurrences where misspecification is stongly suspected for uniform and time-to-event: 
suspect_u <- c(proto = c("VHX_56_2", "VHX_91_2"), # pv3rs vs proto
               uncom = "VHX_39_2", # uncomputable using the prototype
               pwise = c("VHX_113_6", "VHX_450_8", "VHX_489_4", "VHX_529_4", "VHX_532_4")) 
suspect_t <- c(proto = c("VHX_56_2")) 

# Identify outliers using visual inspection:
outlier_u_log <- rhats_tot > 0.42 & (1-Uniform_Pv3Rs[names(rhats_tot),"I"]) < 0.7

if (Figs) pdf(file = "../Figures/compare_prob_prox.pdf", height = 7, width = 12)
par(mfrow = c(2,1), mar = c(5,5,1,2))
plot(x = rhats_tot, y = 1-Uniform_Pv3Rs[names(rhats_tot),"I"], 
     cex = rhats_tot_n/9, bty = "n", 
     pch = PMQ[names(rhats_tot)] + 16, 
     col = outlier_u_log + 1,
     xlab = "Genetic proximity: total relatedness", 
     ylab = "Relapse plus recrudescence\n probability (uniform prior)")
abline(v = q_tot["95%"], lty = "dashed")
text(x = rhats_tot[suspect_u], y = 1-Uniform_Pv3Rs[suspect_u,"I"], 
     labels = suspect_u, pos = 4, cex = 0.5)
legend("right", legend = c("no PQ", "PQ +"), pch = c(16, 17), 
       bty = "n", inset = 0.25)

plot(x = rhats_tot, y = 1-TimeToEvent_Pv3Rs[names(rhats_tot),"I"], 
     cex = rhats_tot_n/9, bty = "n", 
     pch = PMQ[names(rhats_tot)] + 16, 
     col = outlier_u_log + 1,
     xlab = "Genetic proximity: total relatedness", 
     ylab = "Relapse plus recrudescence\n probability (informative prior)")
abline(v = q_tot["95%"], lty = "dashed")
text(x = rhats_tot[suspect_t], y = 1-TimeToEvent_Pv3Rs[suspect_t,"I"], 
     labels = suspect_t, pos = 4, cex = 0.5)
legend("right", legend = c("no PQ", "PQ +"), pch = c(16, 17), 
       bty = "n", inset = 0.25)
if (Figs) dev.off()

if (Figs) pdf(file = "../Figures/outlier_correction.pdf", height = 3.5, width = 12)
par(mfrow = c(1,1), mar = c(5,5,1,2))
X <- Uniform_Pv3Rs[,c("C","L","I")]
X[suspect_u[grepl("pwise", names(suspect_u))],] <-
  Uniform_pwise[suspect_u[grepl("pwise", names(suspect_u))], c("C","L","I")]
plot(x = rhats_tot, y = 1-X[names(rhats_tot),"I"], 
     cex = rhats_tot_n/9, bty = "n", 
     pch = PMQ[names(rhats_tot)] + 16,
     col = c("lightgrey","black")[PMQ[names(rhats_tot)]+1],
     xlab = "Genetic proximity: total relatedness", 
     ylab = "Relapse plus recrudescence\n probability (uniform prior)")
abline(v = q_tot["95%"], lty = "dashed")
outlier_PMQ <- names(which(rhats_tot > 0.45 & (1-X[names(rhats_tot),"I"]) < 0.1 & 
                             PMQ[names(rhats_tot)]))
text(x = rhats_tot[outlier_PMQ], y = 1-Uniform_Pv3Rs[outlier_PMQ,"I"], 
     labels = outlier_PMQ, pos = 4, cex = 0.5)
arrows(x0 = rhats_tot[outlier_PMQ], x1 = rhats_tot[outlier_PMQ], 
       y0 = 1-Uniform_Pv3Rs[outlier_PMQ,"I"] + 0.05, 
       y1 = 0.95, length = 0.05)
legend("right", legend = c("no PQ", "PQ +"), col = c("lightgrey","black"),
       pch = c(16, 17), 
       bty = "n", inset = 0.25)
if (Figs) dev.off()

plot(x = rhats_tot, y = 1-TimeToEvent_Pv3Rs[names(rhats_tot),"I"], 
     cex = rhats_tot_n/9, bty = "n", 
     pch = PMQ[names(rhats_tot)] + 16, 
     col = c("lightgrey","black")[PMQ[names(rhats_tot)]+1],
     xlab = "Genetic proximity: total relatedness", 
     ylab = "Relapse plus recrudescence\n probability (informative prior)")
abline(v = q_tot["95%"], lty = "dashed")
text(x = rhats_tot[outlier_PMQ], y = 1-TimeToEvent_Pv3Rs[outlier_PMQ,"I"], 
     labels = outlier_PMQ, pos = 4, cex = 0.5)
arrows(x0 = rhats_tot[outlier_PMQ], x1 = rhats_tot[outlier_PMQ], 
       y0 = 1-TimeToEvent_Pv3Rs[outlier_PMQ,"I"] + 0.05, 
       y1 = 0.95, length = 0.05)
legend("right", legend = c("no PQ", "PQ +"), col = c("lightgrey","black"),
       pch = c(16, 17), 
       bty = "n", inset = 0.25)


#===============================================================================
# Numbers for ms: 
#===============================================================================
# Average probability of relapse among reinfection-classified 
sum(Uniform_Pv3Rs[which(rhats_tot <= q_tot["95%"]), "L"])/sum(rhats_tot <= q_tot["95%"])
sum(TimeToEvent_Pv3Rs[which(rhats_tot <= q_tot["95%"]), "L"])/sum(rhats_tot <= q_tot["95%"])

# Aside: average probability of recrudescence among reinfection-classified 
sum(Uniform_Pv3Rs[which(rhats_tot <= q_tot["95%"]), "C"])/sum(rhats_one <= q_one["95%"])
sum(TimeToEvent_Pv3Rs[which(rhats_tot <= q_tot["95%"]), "C"])/sum(rhats_one <= q_one["95%"])

# Average probability of reinfection among relapse/recrudescence classified
sum(Uniform_Pv3Rs[which(rhats_tot > q_tot["95%"]), "I"])/sum(rhats_tot > q_tot["95%"])
sum(X[which(rhats_tot > q_tot["95%"]), "I"])/sum(rhats_tot > q_tot["95%"]) # joint-to-pwise outlier correction
sum(TimeToEvent_Pv3Rs[which(rhats_tot > q_tot["95%"]), "I"])/sum(rhats_tot > q_tot["95%"])


#===============================================================================
# Inspect data on major outliers: 
#===============================================================================
outlier_u <- names(which(outlier_u_log))
outlier_extra <- outlier_u[!outlier_u %in% suspect_u]
pids <- unique(apply(do.call(rbind, strsplit(outlier_extra, split = "_"))[,1:2], 1, paste, collapse = "_")) 

# Extract all episodes for pids 
episodes <- unname(unlist(sapply(pids, function(pid) {
  sapply(names(ys_VHX_BPD[[pid]]), function(epi) {
    paste(pid, epi, sep = "_")}, simplify = F)}, simplify = F)))

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
if(Figs) png("../Figures/data_outliers.png", width = 10, height = 7, units = "in", res = 300) 
par(fig = c(0,1,0.2+0.01,1), mar = main_mar, mfrow = c(1,1)) # Important to call before and after plot_data
plot_data(ys = ys_VHX_BPD[pids], fs = fs_VHX_BPD, marker.annotate = F, mar = main_mar, person.vert = T)
par(fig = c(0,1,0.2+0.01,1), mar = main_mar) # Important to call before and after plot_data
text(y = rep(0.0375, length(episodes)),
     x = seq(z, 1-z, length.out = length(episodes)),
     labels = text, cex = 0.4, col = text_col, srt = 90)
points(y = rep(0.06, length(unlist(outlier_extra))),
       x = seq(z, 1-z, length.out = length(episodes))[episodes %in% unlist(outlier_extra)],
       pch = 25, cex = 0.4, bg = "black")
if(Figs) dev.off()


