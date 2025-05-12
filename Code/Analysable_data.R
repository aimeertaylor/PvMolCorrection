################################################################################
# Generate a table of MOIs: sufficient statistic for feasibility of estimation
################################################################################
list(list = ls())
library(Pv3Rs)

# Get vectors of MOIs for all study participants
MOIs <- sapply(ys_VHX_BPD, function(y){
  determine_MOIs(y, return.names = T)
})

# Get maximum number of episodes
n_epi <- sapply(MOIs, function(x) max(as.numeric(names(x))))
max_n_epi <- max(n_epi)
z <- rep(NA, max_n_epi)
names(z) <- paste("Episode", 1:max_n_epi, sep = " ")
tab_MOIs <- t(sapply(MOIs, function(x) {
  z[as.numeric(names(x))] <- x
  return(z)}))

# Add cumulative MOI count and episode count
tab_MOIs <- cbind(tab_MOIs, Total = rowSums(tab_MOIs, na.rm = T), "Episode count" = n_epi)

# For overleaf, compute number of recurrences and participants with paired data
# (length(MOIs) accounts for n_epi = 1 cases because these feature only once in
# length(unlist(MOIs)))
n_recur <- length(unlist(MOIs)) - length(MOIs)
n_pids <- sum(tab_MOIs[,"Total"] > 1)

# For overleaf, compute number of participants with data that can be analysed jointly
n_pids_joint_proto <- sum(tab_MOIs[,"Total"] > 1 & tab_MOIs[,"Total"] < 6 & tab_MOIs[,"Episode count"] < 4)
n_pids_joint_pv3Rs <- sum(tab_MOIs[,"Total"] > 1 & tab_MOIs[,"Total"] < 9)
(n_pids_joint_pv3Rs - n_pids_joint_proto)/n_pids_joint_proto

# Generate excel file to generate figure
write.csv2(tab_MOIs[tab_MOIs[,"Total"] > 5 | tab_MOIs[,"Episode count"] > 3, ], file = "tab_MOIs.csv")

# In excel:
# order rows by MOI total and n_epi
# colour recurrences for which the prototype was not able to generate estimates:
# "VHX_239_2" "VHX_461_2" "VHX_39_2" "VHX_52_2" "VHX_583_2" "VHX_33_2"
# colour participants for which estimates had to be generated pairwise: > 3
# episodes for prototype, > 8 MOI total for Pv3Rs
# save as excel and take screen shot
