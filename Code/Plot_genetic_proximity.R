rm(list = ls())
library(Pv3Rs)
load("../RData/genetic_proximities.RData")
load("../RData/quantiles_null.RData")
load("../RData/marg_results_Pv3Rs.RData") # Load posterior probabilities
Figs <- T

#===============================================================================
# Plot probabilities against genetic proximity: 
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
# Compare with prototype. 
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
# Pv3Rs possible worse for BPD_70_2, VHX_298_2, BPD_253_3  
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
# Plot for ms
par(mfrow = c(2,1), mar = c(5,5,1,2))
#===============================================================================
outlier_u <- c(proto = c("VHX_56_2", "VHX_91_2"), # pv3rs vs proto
               uncom = "VHX_39_2", # uncomputable using the prototype
               pwise = c("VHX_113_6", "VHX_450_8", "VHX_489_4", "VHX_529_4", "VHX_532_4")) 

outlier_t <- c(proto = c("VHX_56_2", "BPD_253_3"), pwise = c("VHX_419_6")) 

plot(x = rhats_tot, y = 1-Uniform_Pv3Rs[names(rhats_tot),"I"], 
     cex = rhats_tot_n/9, bty = "n", 
     pch = PMQ[names(rhats_tot)] + 16, 
     xlab = "Genetic proximity: total relatedness", 
     ylab = "Relapse plus recrudescence\n probability (uniform prior)")
abline(v = q_tot["95%"], lty = "dashed", col = "red")
text(x = rhats_tot[outlier_u], y = 1-Uniform_Pv3Rs[outlier_u,"I"], 
     labels = outlier_u, pos = 4, cex = 0.5)

plot(x = rhats_tot, y = 1-TimeToEvent_Pv3Rs[names(rhats_tot),"I"], 
     cex = rhats_tot_n/9, bty = "n", 
     pch = PMQ[names(rhats_tot)] + 16, 
     xlab = "Genetic proximity: total relatedness", 
     ylab = "Relapse plus recrudescence\n probability (informative prior)")
abline(v = q_tot["95%"], lty = "dashed", col = "red")
text(x = rhats_tot[outlier_t], y = 1-TimeToEvent_Pv3Rs[outlier_t,"I"], 
     labels = outlier_t, pos = 3, cex = 0.5)

plot(x = rhats_tot[names(PMQ[PMQ])], y = 1-Uniform_Pv3Rs[names(PMQ[PMQ]),"I"], 
     cex = rhats_tot_n[names(PMQ[PMQ])]/9, bty = "n", pch = 17, 
     xlim = range(rhats_tot),
     xlab = "Genetic proximity: total relatedness", 
     ylab = "Relapse plus recrudescence\n probability (uniform prior)")
abline(v = q_tot["95%"], lty = "dashed", col = "red")
text(x = rhats_tot[outlier_t[outlier_t %in% names(PMQ[PMQ])]], 
     y = 1-Uniform_Pv3Rs[outlier_t[outlier_t %in% names(PMQ[PMQ])],"I"], 
     labels = outlier_t[outlier_t %in% names(PMQ[PMQ])], pos = 4, cex = 0.5)
outlier <- names(which(rhats_tot[names(PMQ[PMQ])] > 0.4 & 
                         (1-Uniform_Pv3Rs[names(PMQ[PMQ]), "I"]) < 0.8))
text(x = rhats_tot[outlier], y = 1-Uniform_Pv3Rs[outlier,"I"], 
     labels = outlier, pos = 4, cex = 0.5)


plot(x = rhats_tot[names(PMQ[PMQ])], y = 1-TimeToEvent_Pv3Rs[names(PMQ[PMQ]),"I"], 
     cex = rhats_tot_n[names(PMQ[PMQ])]/9, bty = "n", pch = 17, 
     xlim = range(rhats_tot),
     xlab = "Genetic proximity: total relatedness", 
     ylab = "Relapse plus recrudescence\n probability (uniform prior)")
abline(v = q_tot["95%"], lty = "dashed", col = "red")
text(x = rhats_tot[outlier_t[outlier_t %in% names(PMQ[PMQ])]], 
     y = 1-TimeToEvent_Pv3Rs[outlier_t[outlier_t %in% names(PMQ[PMQ])],"I"], 
     labels = outlier_t[outlier_t %in% names(PMQ[PMQ])], pos = 4, cex = 0.5)
text(x = rhats_tot[outlier], y = 1-TimeToEvent_Pv3Rs[outlier,"I"], 
     labels = outlier, pos = 4, cex = 0.5)


#===============================================================================
# Numbers of ms: 
#===============================================================================
# Average probability of relapse among reinfection-classified 
sum(Uniform_Pv3Rs[which(rhats_tot <= q_tot), "L"])/sum(rhats_one <= q_one)
sum(TimeToEvent_Pv3Rs[which(rhats_tot <= q_tot), "L"])/sum(rhats_one <= q_one)

# Aside: average probability of recrudescence among reinfection-classified 
sum(Uniform_Pv3Rs[which(rhats_tot <= q_tot), "C"])/sum(rhats_one <= q_one)
sum(TimeToEvent_Pv3Rs[which(rhats_tot <= q_tot), "C"])/sum(rhats_one <= q_one)

# Proportion / probability of "relapse or recrudescence" recurrences among
# recurrences in PMQ treated participants:
sum((1-Uniform_Pv3Rs[,"I"])*PMQ[rownames(Uniform_Pv3Rs)])/sum(PMQ[rownames(Uniform_Pv3Rs)])
sum((1-TimeToEvent_Pv3Rs[,"I"])*PMQ[rownames(TimeToEvent_Pv3Rs)])/sum(PMQ[rownames(TimeToEvent_Pv3Rs)])
sum((rhats_tot > q_tot)*PMQ[names(rhats_tot)])/sum(PMQ[names(rhats_tot)])

# Proportion / probability of "relapse or recrudescence" recurrences among
# reinfection-classified recurrences in PMQ treated participants:
sum((1-Uniform_Pv3Rs[,"I"])*PMQ[pids]*(rhats_one <= q_one))/sum(PMQ[pids]*(rhats_one <= q_one))
sum((1-TimeToEvent_Pv3Rs[,"I"])*PMQ[pids]*(rhats_one <= q_one))/sum(PMQ[pids])

#===============================================================================
# Inspect data on outliers: 
#===============================================================================
pids <- unique(apply(do.call(rbind, strsplit(FN_tot, split = "_"))[,1:2], 1, paste, collapse = "_")) 

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


