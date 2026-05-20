################################################################################
# 
################################################################################
library(Pv3Rs)
source("compute_posterior_approxjoint.R")

# Synthetic data: (check Understanding posterior)
yhet <- list("0" = list(m1 = "1", m2 = "1", m3 = "1"),
          "1" = list(m1 = "2", m2 = "2", m3 = "2"),
          "2" = list(m1 = "3", m2 = "3", m3 = "3"),
          "3" = list(m1 = "4", m2 = "4", m3 = "4"),
          "4" = list(m1 = "5", m2 = "5", m3 = "5"))
# Synthetic data: (check Understanding posterior)
yhom <- list("0" = list(m1 = "1", m2 = "1", m3 = "1"),
             "1" = list(m1 = "1", m2 = "1", m3 = "1"),
             "2" = list(m1 = "1", m2 = "1", m3 = "1"),
             "3" = list(m1 = "1", m2 = "1", m3 = "1"),
             "4" = list(m1 = "1", m2 = "1", m3 = "1"))
fs <- list(m1 = c('1' = 0.2, '2' = 0.2, '3' = 0.2, '4' = 0.2, '5' = 0.2),
           m2 = c('1' = 0.2, '2' = 0.2, '3' = 0.2, '4' = 0.2, '5' = 0.2),
           m3 = c('1' = 0.2, '2' = 0.2, '3' = 0.2, '4' = 0.2, '5' = 0.2))
suppressMessages(compute_posterior(yhet, fs, progress.bar = FALSE))$marg
suppressMessages(compute_posterior_approxjoint(yhet, fs))
suppressMessages(compute_posterior(yhom, fs, progress.bar = FALSE))$marg
suppressMessages(compute_posterior_approxjoint(yhom, fs))

# Look at differences between joint and approximate joint:
x <- sort(abs(marg_pwise[["Uniform"]][, "I"] - marg_joint[["Uniform"]][, "I"]))
y <- as.numeric(do.call(rbind, strsplit(names(x), split = "_"))[,3])
plot(x = y, y = x)

# Choose pid and plot data:
#"VHX_66" - not supporting all sibling graph
#"VHX_621" (some matches and a relapse) and "BPD_455" (some matches) and "BPD_234" (no matches) - increasing graph size
pid <- "VHX_621"
plot_data(ys_VHX_BPD[pid], fs = sapply(fs_VHX_BPD, sort, decreasing = T), 
          marker.annotate = F)

# Compute recurrence state probabilities:
joint <- compute_posterior(ys_VHX_BPD[[pid]], fs = fs_VHX_BPD, return.logp = T)          
aprox <- compute_posterior_approxjoint(ys_VHX_BPD[[pid]], fs = fs_VHX_BPD)

# Inspect recurrence state probabilities:
joint$marg
aprox
suppressMessages(compute_posterior(ys_VHX_BPD[[pid]][c(1,2)], fs = fs_VHX_BPD))$marg
suppressMessages(compute_posterior(ys_VHX_BPD[[pid]][c(1,3)], fs = fs_VHX_BPD))$marg
suppressMessages(compute_posterior(ys_VHX_BPD[[pid]][c(2,3)], fs = fs_VHX_BPD))$marg

RG_liks <- sapply(joint$RGs, function(x) x$logp)
RG_inds <- which(!is.infinite(RG_lik))
mois <- determine_MOIs(ys_VHX_BPD[[pid]])
par(mfrow = c(3,2))
sapply((1:length(RG_liks)), function(ind) {
  RG <- joint$RGs[[ind]]
  plot_RG(RG_to_igraph(RG, mois), edge.curved = 0.5)
  title(RG_liks[ind])
  box()})

# For VHX_66, I think if the likelihood of the all sibling graph were not -Inf,
# the probabilities of IL and II would be equal whereas:
sort(joint$joint, decreasing = T)
abs(joint$marg[,"I"] - aprox[,"I"])
      