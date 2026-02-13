################################################################################
# Functions required to compute genetic proximity with dcifer
#
# From vignetteDcifer: If afreq contains “extra” alleles that are not listed in
# dsmp, these alleles are added to dsmp. If dsmp has alleles not included in
# afreq or listed with a frequency of 0, they are assigned a small probability
# (optional minfreq argument) with subsequent renormalization.
################################################################################

# List marker and per-marker allele names 
markers <- names(fs_VHX_BPD)
alleles_per_marker <- sapply(fs_VHX_BPD, names)

# Function to convert per-episode data formatted for pv3Rs (categorical lists
# per marker) into dcifer format (binary list across all alleles per marker)
pv3rs_to_dcifer <- function(y_epi, alleles_per_marker){
  z <- sapply(markers, function(m) {
    x <- rep(0, length(alleles_per_marker[[m]])) # construct binary vector
    names(x) <- alleles_per_marker[[m]] # name binary vector with alleles
    if (m %in% names(y_epi)) { # If m listed for y_epi
      alleles <- y_epi[[m]] # extract alleles (NAs possibly included)
      if (!all(is.na(alleles))) { # If m listed and successfully typed
        x[as.character(alleles[!is.na(alleles)])] <- 1 # add one for non-na detected alleles
      }}
    return(x)
  })
  return(z)
}

# Function returns for each recurrence maximum genetic proximity using dcifer
# relative to all previous episodes per participants
genetic_proximity <- function(y, type = "rhat_one"){ # per participant
  
  n_epi <- length(y) # if at least one recurrence
  
  if(n_epi > 1) { 
    max_rhat <- sapply(2:n_epi, function(i){ # max_rhat with max_nmarks and episode name
      prev_rhats_nmarks <- sapply(1:(i-1), function(j) { # compare to all previous
        # extract episode name (some ys have data on non-sequential episodes, e.g., VHX_419)
        episode <- names(y)[j]
        # extract paired episode data in dcifer format: 
        pair <- sapply(y[c(i,j)], pv3rs_to_dcifer, alleles_per_marker, simplify = F) 
        # compute number of markers with comparible data for pair ij
        nmarks <- sum(apply(sapply(pair, function(epi) sapply(epi, sum)) > 0, 1, all))
        
        if(nmarks == 0) {
          rhat = NA
        } else {
          cois <- determine_MOIs(y[c(i,j)]) # Use same MOIs as Pv3Rs 
          afreq <- fs_VHX_BPD # Use same allele frequencies as Pv3Rs
          if (type == "rhat_one") { # Assuming at most one genotype related 
            rhat <- dcifer::ibdPair(pair, cois, afreq, M = 1)            
          } else if (type == "rhat_avg") { # Quick approximate average relatedness
            rhat <- unique(dcifer::ibdEstM(pair, cois, afreq, equalr = TRUE)) 
          } else if (type == "rhat_tot") {# Quick approximate total relatedness
            rhat <- sum(dcifer::ibdEstM(pair, cois, afreq, equalr = TRUE)) 
          }
        }
        return(list(prev_epi = episode, rhat = rhat, nmarks = nmarks))
      })
      
      # Extract maximum rhat; don't use which.max since possibly multiple
      # comparisons with max value in which case we want the one with largest
      # nmarks, not first (the return value which.max)
      rhat_max <- max(unlist(prev_rhats_nmarks["rhat", ]), na.rm = T) 
      rhat_max_i <- which(unlist(prev_rhats_nmarks["rhat", ]) == rhat_max)
      
      if (length(rhat_max_i) > 1) { # If multiple rhat_max values
        # If among multiple rhat_max equally high nmarks, return first using which.max: 
        nmarks_max_i <- which.max(unlist(prev_rhats_nmarks["nmarks", rhat_max_i])) 
        return(prev_rhats_nmarks[, rhat_max_i[nmarks_max_i]])
      } else {
        return(prev_rhats_nmarks[, rhat_max_i])
      }
      
    })
    colnames(max_rhat) <- names(y)[2:n_epi]
    return(max_rhat)
  } else {
    return(NA)
  }
}

#===============================================================================
# Sandbox: testing some dcifer calculations on VHX BPD
#===============================================================================
# pid <- "VHX_532"
# y <- ys_VHX_BPD[[pid]]
# plot_data(ys_VHX_BPD[pid], fs = fs_VHX_BPD, marker.annotate = F)
# i <- 5; j = 1
# pair <- sapply(y[c(i,j)], pv3rs_to_dcifer, alleles_per_marker, simplify = F)
# cois <- determine_MOIs(y[c(i,j)])
# revals <- mapply(generateReval, 1:5, nr = c(1e3, 1e2, 32, 16, 12))
# # Compute r values for each genotype pair up to max(cois) assuming pairs are differentially related
# rdiff <- ibdEstM(pair, cois, afreq = fs_VHX_BPD, Mmax = max(cois), equalr = FALSE, reval = revals) # slow
# # Compute r values for each genotype pair up to max(cois) assuming pairs are equally related
# rsame <- ibdEstM(pair, cois, afreq = fs_VHX_BPD, equalr = TRUE) # quick
# # Average of different r values is close to the repeated value when genotypes are considered equal: 
# mean(rdiff); rsame[1] 

#===============================================================================
# Close correlation between average relatedness and approximate average
# relatedness is illustrated in code copied from vignetteDcifer:
#===============================================================================
# afile <- system.file("extdata", "MozAfreq.csv", package = "dcifer")
# afreq <- readAfreq(afile, lvar = "locus", avar = "allele", fvar = "freq")
# meta <- unique(read.csv(sfile)[c("sampleID", "province")])
# meta <- meta[match(names(dsmp), meta$sampleID), ] # order samples as in dsmp
# provinces <- c("Maputo", "Inhambane")
# nsite <- table(meta$province)[provinces]
# ord <- order(factor(meta$province, levels = provinces))
# dsmp <- dsmp[ord]
# coi <- coi[ord]
# dres <- ibdDat(dsmp, coi, afreq, pval = TRUE, confint = TRUE, rnull = 0, alpha = 0.05, nr = 1e3) # Assumes M = 1
# alpha <- 0.05 # significance level
# sig <- dres[, , "p_value"] <= alpha
# isig <- which(sig, arr.ind = TRUE)
# nsig <- nrow(isig)
# sig1 <- sig2 <- vector("list", nsig)
# for (i in 1:nsig) {
#   sig1[[i]] <- ibdEstM(dsmp[isig[i, ]], coi[isig[i, ]], afreq, Mmax = 5,
#                        equalr = FALSE, reval = revals) 
# }
# for (i in 1:nsig) {
#   sig2[[i]] <- ibdEstM(dsmp[isig[i, ]], coi[isig[i, ]], afreq, equalr = TRUE)
# }
# M1 <- sapply(sig1, function(r) sum(r > 0))
# M2 <- sapply(sig2, length)
# rtotal1 <- sapply(sig1, sum)
# rtotal2 <- sapply(sig2, sum)

