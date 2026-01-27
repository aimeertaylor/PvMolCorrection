################################################################################
# Generate a table that contains per-episode MOIs: sufficient statistic for
# feasibility of estimation using both the prototype and Pv3Rs
################################################################################
list(list = ls())
library(Pv3Rs) # For plots
Figs <- FALSE
  
# Load data from all participants to extract meta data
load("../jwatowatson-RecurrentVivax-4870715/RData/TimingModel/Combined_Time_Event.RData") 

# Episode counts and treatment for all participant IDs inc. those whose episodes were not typed
epi_counts <- sapply(split(Combined_Time_Data$episode[Combined_Time_Data$Censored != 1], 
                    Combined_Time_Data$patientid[Combined_Time_Data$Censored != 1]), max)
arm <- sapply(split(Combined_Time_Data$arm_num, Combined_Time_Data$patientid), unique) 
pmq <- sapply(arm, function(x) grepl("PMQ", x))

# Get vectors of MOIs for all study participants with one ore more typed episode
MOIs <- sapply(Pv3Rs::ys_VHX_BPD, function(y) determine_MOIs(y, return.names = T))

# Get MOI summaries for typed episodes
typd_epi_count <- sapply(MOIs, function(x) length(as.numeric(names(x))))
max_typd_epi_count <- max(typd_epi_count) # Maximum episode count typed 
z <- rep(NA, max_typd_epi_count); names(z) <- paste("Episode", 1:max_typd_epi_count, sep = " ")
tab_MOIs <- t(sapply(MOIs, function(x) {
  z[as.numeric(names(x))] <- x
  return(z)}))

# Add cumulative MOI count  count
tab_MOIs <- data.frame(tab_MOIs, Total = rowSums(tab_MOIs, na.rm = T), 
                       "Typed.episode.count" = typd_epi_count, 
                       "Total.episode.count" = epi_counts[names(typd_epi_count)],
                       "Treatment.received" = arm[names(typd_epi_count)])


# For overleaf, compute number of recurrences and participants with paired data
# (length(MOIs) accounts for typd_epi_count = 1 cases because these feature only once in
# length(unlist(MOIs)))
n_recur <- length(unlist(MOIs)) - length(MOIs) # No. of recurrences with paired data
n_pids <- sum(tab_MOIs$Typed.episode.count > 1) # No. of pids with paired data

# For overleaf, compute number of participants with data that can be analysed jointly
n_pids_joint_proto <- sum(tab_MOIs$Typed.episode.count > 1 & tab_MOIs$Typed.episode.count < 4 & tab_MOIs$Total < 6)
n_pids_joint_pv3Rs <- sum(tab_MOIs$Typed.episode.count > 1 & tab_MOIs$Total < 9)

# For overleaf, compute the % decrease in number of pids whose data were modelled pairwise 
((n_pids-n_pids_joint_pv3Rs) - (n_pids-n_pids_joint_proto-4))/(n_pids-n_pids_joint_proto-4)

# Generate excel file to generate figure
write.csv2(tab_MOIs[tab_MOIs$Total > 5 | tab_MOIs$Typed.episode.count > 3, ], 
           file = "../Figures/tab_MOIs.csv")

# In excel:
# order rows by MOI total, epi_counts, and then typd_epi_count
# colour recurrences for which the prototype was not able to generate estimates:
# "VHX_239_2" "VHX_461_2" "VHX_39_2" "VHX_52_2" "VHX_583_2" "VHX_33_2"
# colour participants for which estimates had to be generated pairwise: 
# > 3 episodes for prototype, 
# > 8 MOI total for Pv3Rs,
# save as excel and take screen shot

# Plot all paired data
if(Figs) png("../Figures/data_paired.png", width = 7, height = 7, units = "in", res = 300)
Pv3Rs::plot_data(ys = ys_VHX_BPD[names(which(typd_epi_count > 1))], 
                 fs = fs_VHX_BPD, marker.annotate = F, person.vert = T,
                 mar = c(2,2.8,1,2.8)) # Margin around main plot)
if(Figs) dev.off()

# Plot data for pids with a recurrence whose recurrent state probability could
# not be estimated using the prototype; annotate with probabilities given uniform prior
pids <- c("VHX_239","VHX_461","VHX_39","VHX_52","VHX_583","VHX_33")
load("../RData/results_Pv3Rs.RData")
marg <- sapply(pids, function(pid) {
  if(is.na(Uniform_joint[pid])) {
    Uniform_pairwise[[pid]]
  } else {
    Uniform_joint[[pid]][["marg"]]
  }
})
state_names <- c("Recrud.", "Relapse", "Reinf.")
episodes <- unlist(sapply(marg, function(x) c(1, as.numeric(rownames(x)))))
probs <- unlist(sapply(marg, function(x) apply(x, 1, max))) # Get max marginal estimates 
states <- unlist(sapply(marg, function(x) apply(x, 1, function(z) state_names [which.max(z)])))
main_mar <- c(5,2.6,1,2.6) # Margin around main plot


if(Figs) png("../Figures/data_unestimatable.png", width = 7, height = 7, units = "in", res = 300)
Pv3Rs::plot_data(ys = Pv3Rs::ys_VHX_BPD[pids], mar = main_mar, 
                 fs = fs_VHX_BPD, marker.annotate = F)
par(fig = c(0,1,0,1), mar = main_mar) # Reset before text annotation (important)
#box() # Toggle on and off to see plot limits for text() — not enacted in margin
text(y = rep(0.21, length(probs)), adj = 1, # right justify 
     x = seq(0.025, 0.98, length.out = length(episodes))[-which(episodes == 1)], 
     labels = paste0(states, " ", 100*round(probs,2), "%"), 
     cex = 0.7, srt = 90)
if(Figs) dev.off()


