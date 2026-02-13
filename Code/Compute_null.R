################################################################################
# To do: back compute alpha that would give same failure rate as that previously 
# published
################################################################################

rm(list = ls())
load("../RData/ys_null.RData")
source("genetic_proximity.R")
source("match_counting.R")
load("../RData/probs_null.RData") # Pv3Rs probabilities for null data
alphas <- c(0.05) # Set equal to the FDR of Pv3Rs
run_null <- F

# For each pair compute number of markers with allele data on both episodes:
n_marks <- sapply(ys_null, function(y){
  n_mark <- sum(sapply(y[[1]], function(x) !all(is.na(x))) & 
                  sapply(y[[2]], function(x) !all(is.na(x))))  
})
one_plus <- names(which(n_marks > 0)) # not all comparisons have comparable data
n_one_plus <- length(one_plus) # 247873 comparisons with data one plus marker

if(run_null){
  tictoc::tic()
  null_one <- sapply(ys_null[one_plus], function(y) unlist(genetic_proximity(y, type = "rhat_one")["rhat",]))
  null_avg <- sapply(ys_null[one_plus], function(y) unlist(genetic_proximity(y, type = "rhat_avg")["rhat",]))
  null_tot <- sapply(ys_null[one_plus], function(y) unlist(genetic_proximity(y, type = "rhat_tot")["rhat",]))
  null_mat <- sapply(ys_null[one_plus], function(y) match_counting(y)["p",])
  tictoc::toc()
  save(null_one, null_avg, null_mat, null_tot, file = "../RData/nulls.RData")  
} else {
  load("../RData/nulls.RData")
}

plot(x = null_one, y = null_avg, bty = "n", pch = 20)
abline(a = 0, b = 1)
q_one <- quantile(null_one, probs = 1-alphas)
q_avg <- quantile(null_avg, probs = 1-alphas)
q_tot <- quantile(null_avg, probs = 1-alphas)
q_mat <- quantile(null_avg, probs = 1-alphas)
q_3Rs <- quantile(x = 1-probs[,"I"], probs = 1-alphas)
par(mfrow = c(2,2))
hist(null_one); abline(v = q_one, col = "red")
hist(null_avg); abline(v = q_avg, col = "red")
hist(null_tot); abline(v = q_tot, col = "red")
hist(null_mat); abline(v = q_mat, col = "red")
hist(1-probs[,"I"]); abline(v = q_3Rs, col = "red")

save(q_one, q_avg, q_mat, q_3Rs, q_tot, file = "../RData/quantiles_null.RData")
