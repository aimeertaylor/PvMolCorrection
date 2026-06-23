################################################################################
# Purpose: compute the probability of allele matches under the assumption of iid
# draws parameterised by allele frequencies and compare with Dcifer and Pv3Rs
#
# Results: strong trend with no major outliers between probabilities of chance
# matching and Dcifer. As such, comparison of chance matching and Pv3Rs
# generates same insight as Dcifer vs Pv3Rs. Chance matching is easier to inpret
# than Dcifer but Dcifer is the current state-of-the-art. As such, could be
# viewed as more valuable for the community to compare with Dcifer.
################################################################################
library(Pv3Rs)
rm(list = ls())

#===============================================================================
# Function to compute the probability of the observed matches assuming all
# alleles are independent and identically draws from per-marker distributions
#===============================================================================
chance_matching <- function(x){
  n_epi <- length(x)
  if(n_epi > 1) {
    
    markers <- names(x[[1]])
    n_markers <- length(markers)
    
    if (!all(sapply(x, function(x_epi) all(names(x_epi) == markers)))) {
      stop("Not all episodes typed at the same markers")}
    
    z <- sapply(2:n_epi, function(i){ # For each recurrence
      
      z_all <- sapply(1:(i-1), function(j) { # Compare to all previous episodes
        
        # Get matching alleles:
        mat_all <- sapply(markers, function(m) intersect(x[[i]][[m]], x[[j]][[m]]), simplify = F)
        
        if(sum(unlist(mat_all), na.rm = T) == 0){ # If no matching alleles
          iid <- 1
        } else {
          mat_frq <- unlist(sapply(markers, function(m) fs_VHX_BPD[[m]][as.character(mat_all[[m]])]))
          iid <- sum(log(mat_frq^2), na.rm = T)
        }
        return(iid)
      })
      return(min(z_all, na.rm = T))
    })
    
    names(z) <- names(x)[-1] # Return episode names if available
    return(z)
  } else {
    stop("Unpaired data")
  }
}

#===============================================================================
# Some sanity checks
#===============================================================================
fs_VHX_BPD[["PV.3.502"]]
fs_VHX_BPD[["PV.ms8"]]

x1 <- list(initial = list("PV.3.502" = NA, "PV.ms8" = NA), 
           recur = list("PV.3.502" = NA, "PV.ms8" = NA), 
           recur2 = list("PV.3.502" = NA, "PV.ms8" = NA))

x2 <- list(initial = list("PV.3.502" = c(1,2), "PV.ms8" = c(1,4)), 
           recur = list("PV.3.502" = 1, "PV.ms8" = 1))

x3 <- list(initial = list("PV.3.502" = c(3,2), "PV.ms8" = c(1,4)), 
           recur = list("PV.3.502" = 1, "PV.ms8" = 1))

x4 <- list(initial = list("PV.3.502" = c(3,1), "PV.ms8" = 3), 
           recur = list("PV.3.502" = 2, "PV.ms8" = c(1,2)))

chance_matching(x1)
chance_matching(x2)
chance_matching(x3)
chance_matching(x4)

#===============================================================================
# Compute matching probabilities for pids with one or more recurrences
#===============================================================================
n_epi_per_pid <- sapply(ys_VHX_BPD, length)
log_prob_iid <- unlist(sapply(ys_VHX_BPD[n_epi_per_pid > 1], chance_matching))
names(log_prob_iid) <- gsub(".", "_", names(log_prob_iid), fixed = T)

#===============================================================================
# Compare probabilities assuming IID draws with Dcifer relatedness - agree
# Can toggle rhats_tot with rhats_one and rhats_avg
#===============================================================================
load("../RData/genetic_proximities.RData")
plot(x = log_prob_iid[names(rhats_tot)]/rhats_tot_n, y = rhats_tot, 
     cex = rhats_tot_n/9, pch = 19, bty = "n", 
     xlab = "Per-marker average log (matching probability assuming iid draws)", 
     ylab = "Total relatedness")

#===============================================================================
# Compare probabilities assuming IID draws with posterior probabilities of reinfection
#===============================================================================
load("../RData/marg_results_Pv3Rs.RData") # Load posterior probabilities
plot(x = log_prob_iid[rownames(Uniform_Pv3Rs)]/rhats_tot_n[rownames(Uniform_Pv3Rs)], 
     y = Uniform_Pv3Rs[,"I"], 
     cex = rhats_tot_n[rownames(Uniform_Pv3Rs)]/9, pch = 19, bty = "n", 
     xlab = "Per-marker average log (matching probability assuming iid draws)", 
     ylab = "Reinfection posterior probability")

suspect_u <- c(proto = c("VHX_56_2", "VHX_91_2"), # pv3rs vs proto
               uncom = "VHX_39_2", # uncomputable using the prototype
               pwise = c("VHX_113_6", "VHX_450_8", "VHX_489_4", "VHX_529_4", "VHX_532_4")) 

text(x = log_prob_iid[suspect_u]/rhats_tot_n[suspect_u], y = Uniform_Pv3Rs[suspect_u,"I"], 
     labels = suspect_u, pos = 4, cex = 0.5)


