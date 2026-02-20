################################################################################
# Generate a table that contains per-episode MOIs: sufficient statistic for
# feasibility of estimation using Pv3Rs. It is not a sufficient summary
# statistic for the prototype where limits were not clear cut (Max_Eps = 3 and
# Max_Tot_Vtx = 6 led to attempted joint calculations for four total MOI = 6
# participants: VHX 239, 461 and 52 failed, VHX 214 succeeded; see
# Extract_and_plot_prototype.R)
################################################################################
list(list = ls())
library(Pv3Rs) # For plots
Figs <- T

# Length (in markers) of each episode:
range(unlist(sapply(ys_VHX_BPD, function(y) sapply(y, length))))

# Numbers of typed episodes per person
table(sapply(ys_VHX_BPD, length))[-1]

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
log_joint_proto <- tab_MOIs$Typed.episode.count > 1 & tab_MOIs$Typed.episode.count < 4 & tab_MOIs$Total < 7
names(log_joint_proto) <- rownames(tab_MOIs)
log_joint_proto[c("VHX_239", "VHX_461", "VHX_52")] <- FALSE # failed attempts with total MOI = 6
n_pids_joint_proto <- sum(log_joint_proto)
n_recur_joint_proto <- length(unlist(MOIs[log_joint_proto])) - length(MOIs[log_joint_proto])

log_joint_pv3Rs <- tab_MOIs$Typed.episode.count > 1 & tab_MOIs$Total < 9
n_pids_joint_pv3Rs <- sum(log_joint_pv3Rs) # Includes people with only one recurrence
n_recur_joint_pv3Rs <- length(unlist(MOIs[log_joint_pv3Rs])) - length(MOIs[log_joint_pv3Rs])

# Number of participants where joint inference is pairwise
log_one_recur_proto <- tab_MOIs$Typed.episode.count == 2 & tab_MOIs$Total < 6 # Excludes VHX_214 
log_one_recur_pv3Rs <- tab_MOIs$Typed.episode.count == 2 & tab_MOIs$Total < 9
n_pids_one_recur_proto <- sum(log_one_recur_proto)
n_pids_one_recur_pv3Rs <- sum(log_one_recur_pv3Rs)

# For overleaf (demoted), compute the % decrease in number of pids whose data
# were modelled using an approximation of joint inference
((n_pids-n_pids_joint_pv3Rs) - (n_pids-n_pids_joint_proto-4))/(n_pids-n_pids_joint_proto-4)

# Generate excel file to generate figure
write.csv(tab_MOIs[tab_MOIs$Total > 5 | tab_MOIs$Typed.episode.count > 3, ], 
           file = "../Figures/tab_MOIs.csv")

# In excel:
# order rows by MOI total, epi_counts, and then typd_epi_count
# colour recurrences for which the prototype was not able to generate estimates:
# "VHX_239_2" "VHX_461_2" "VHX_39_2" "VHX_52_2" "VHX_583_2" "VHX_33_2"
# delete row "VHX_214"
# colour participants for which estimates had to be generated pairwise: 
# > 3 episodes for prototype, 
# > 8 MOI total for Pv3Rs,
# save as excel and take screen shot

# Plot all paired data: sort by treatment within trial and in ascending id order
names_repeat_episode <- names(which(typd_epi_count > 1)) # alpha-numeric sorting
names_repeat_episode_BPD <- names_repeat_episode[grepl("BPD", names_repeat_episode)]
names_repeat_episode_VHX <- names_repeat_episode[grepl("VHX", names_repeat_episode)]
names_repeat_episode_VHX_PMQ <- names_repeat_episode_VHX[grepl("PMQ", arm[names_repeat_episode_VHX])]                                  
names_repeat_episode_VHX_noPMQ <- names_repeat_episode_VHX[!grepl("PMQ", arm[names_repeat_episode_VHX])]                                  
names_repeat_episode_arm <- c(names_repeat_episode_BPD, 
                              names_repeat_episode_VHX_PMQ, 
                              names_repeat_episode_VHX_noPMQ)

if(Figs) png("../Figures/data_repeat_episode.png", res = 300, height = 9, width = 20, units = "in")
Pv3Rs::plot_data(ys = ys_VHX_BPD[names_repeat_episode_arm], 
                 fs = fs_VHX_BPD, marker.annotate = F, person.vert = T,
                 mar = c(2,3.5,1,3), 
                 gridlines = F) # Margin around main plot)
if(Figs) dev.off()

# Plot data for pids with a recurrence whose recurrent state probability could
# not be estimated using the prototype; annotate with probabilities given uniform prior
# Order according to order mentioned in overleaf
pids <- c("VHX_461","VHX_239","VHX_33","VHX_52","VHX_583","VHX_39")
load("../RData/results_Pv3Rs.RData")
marg <- sapply(pids, function(pid) {
  if(is.na(Uniform_joint[pid])) {
    Uniform_pairwise[[pid]]
  } else {
    Uniform_joint[[pid]][["marg"]]
  }
})

state_names <- c("C", "L", "I")
episodes <- unname(unlist(sapply(pids, function(pid) {
  sapply(names(ys_VHX_BPD[[pid]]), function(epi) {
    paste(pid, epi, sep = "_")})})))
recurrences <- unlist(sapply(pids, function(pid) paste0(pid, "_", rownames(marg[[pid]]))))
probs <- unlist(sapply(marg, function(x) apply(x, 1, max))) # Get max marginal estimates 
states <- unlist(sapply(marg, function(x) apply(x, 1, function(z) state_names [which.max(z)])))
no_episodes_per_pid <- sapply(marg, nrow) + 1 # Needed for text_col
text <- rep("", length(episodes)); names(text) <- episodes
text[recurrences] <- paste0(states, " ", 100*round(probs,2))
text_col <- unlist(sapply(1:length(pids), function(i){
  if(i %% 2 == 1) {
    rep("white", no_episodes_per_pid[i])
  } else {
    rep("black", no_episodes_per_pid[i])
  }}))
main_mar <- c(1.5, 3.5, 1.5, 3.5)
z <- 1/(2*length(episodes)) # See text x placement

if(Figs) png("../Figures/data_uncomputable.png", res = 300, height = 7, width = 10, units = "in")
par(fig = c(0,1,0.2+0.01,1), mar = main_mar) # Important to call before and after plot_data 
Pv3Rs::plot_data(ys = Pv3Rs::ys_VHX_BPD[pids], 
                fs = fs_VHX_BPD, marker.annotate = F, mar = main_mar)
par(fig = c(0,1,0.2+0.01,1), mar = main_mar) # Important to call before and after plot_data
text(y = rep(0.011, length(episodes)), # right justify 
     x = seq(z, 1-z, length.out = length(episodes)), 
     labels = text, cex = 0.6, col = text_col)
points(y = rep(-0.0075, length(pids)), # Only one recurrence per pid was not computable
       x = seq(z, 1-z, length.out = length(episodes))[episodes %in% paste0(pids, "_2")], 
       pch = 17, cex = 0.5) 
dev.off()

