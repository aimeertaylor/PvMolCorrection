
rm(list = ls())
load("../Rdata/MS_final_generated_by_running_all_chunks_of_Pooled_Analysis.Rmd")

# Create a matrix of estimates 
TimeToEvent_old <- data.frame(C = MS_final[,"C_median"], L = MS_final[,"L_median"], I = MS_final[,"I_median"])
rownames(TimeToEvent_old) <- MS_final[,"Episode_Identifier"]
TimeToEvent_old$joint <- MS_final[,"TotalEpisodes"] <= 3


#===============================================================================
# Compare Time-to-Event results
#===============================================================================
load("../RData/marg_results_Pv3Rs.RData")

# Check overlap 
all(rownames(TimeToEvent_old) %in% rownames(TimeToEvent_Pv3Rs))
all(rownames(TimeToEvent_Pv3Rs) %in% rownames(TimeToEvent_old))

# Check missing meets expectation: yes; see Analysable_data.png
rownames(TimeToEvent_Pv3Rs)[which(!rownames(TimeToEvent_Pv3Rs) %in% rownames(TimeToEvent_old))]
# "VHX_239_2" "VHX_33_2"  "VHX_39_2"  "VHX_461_2" "VHX_52_2"  "VHX_583_2"

# Compare old and new relapse probabilities
par(mfrow = c(2,2))
plot(y = TimeToEvent_old[,"I"], 
     x = TimeToEvent_Pv3Rs[rownames(TimeToEvent_old), "I"],
     pch = 20, col = 1+grepl("BPD", rownames(TimeToEvent_old)))
abline(a = 0, b = 1, lty = "dotted")
plot(y = TimeToEvent_old[,"C"], 
     x = TimeToEvent_Pv3Rs[rownames(TimeToEvent_old), "C"],
     pch = 20, col = 1+grepl("BPD", rownames(TimeToEvent_old)))
abline(a = 0, b = 1, lty = "dotted")

plot(y = TimeToEvent_old[,"L"], 
     x = TimeToEvent_Pv3Rs[rownames(TimeToEvent_old), "L"],
     pch = 20, col = 1+grepl("BPD", rownames(TimeToEvent_old)))
abline(a = 0, b = 1, lty = "dotted")
# Extract and annotate big differences 
diffs <- abs(TimeToEvent_old[,"L"] - TimeToEvent_Pv3Rs[rownames(TimeToEvent_old), "L"])
big_diffs <- rownames(TimeToEvent_old)[which(diffs > 0.4)]
text(x = TimeToEvent_Pv3Rs[big_diffs, "L"], y = TimeToEvent_old[big_diffs, "L"], 
     labels = big_diffs, pos = c(4,2))

# Plot data and inspect estimates for the participants with estimates that differ
# Both 5/9 match; rarer alleles for BPD_562...
plot_data(ys = ys_VHX_BPD[c("VHX_56", "BPD_562")], fs = fs_VHX_BPD)
TimeToEvent_old[big_diffs,] 
TimeToEvent_Pv3Rs[big_diffs,] # BPD more reasonable; VHX: half sib missclassification 

