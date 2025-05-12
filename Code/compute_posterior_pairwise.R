# Pairwise inference 
compute_posterior_pairwise <- function(y, fs, prior = NULL) {

  # Number of episodes 
  n_epi <- length(y)
  
  # Create an array to store results
  z <- array(NA, dim = c(n_epi, n_epi, 3), dimnames = list(names(y), names(y), c("C","L","I")))
  
  # Generate a list of episodes whose data will be analysed pairwise 
  epi_pairs <- gtools::combinations(n = length(names(y)), r = 2, v = as.numeric(names(y)))
  
  for(i in 1:nrow(epi_pairs)) {
    epi_pair_chr <- as.character(epi_pairs[i,])
    y_pair <- y[epi_pair_chr] 
    if(sum(determine_MOIs(y_pair)) < 8) { # Check fewer than 8 genotypes
      if(is.null(prior)) {
        x <- compute_posterior(y = y_pair, fs = fs)$marg             
      } else {
        rec_prior <- prior[epi_pair_chr[2], , drop = F] # Extract prior for single recurrence 
        x <- compute_posterior(y = y_pair, fs = fs, prior = rec_prior)$marg  
      }
 
      z[epi_pair_chr[1],epi_pair_chr[2],"C"] <- x[1,"C"]
      z[epi_pair_chr[1],epi_pair_chr[2],"L"]  <- x[1,"L"]
      z[epi_pair_chr[1],epi_pair_chr[2],"I"]  <- x[1,"I"]
    }
  }
  
  # For reinfection, compare to all preceding episodes
  I_unnorm <- apply(z[,-1,"I", drop = F], 2, min, na.rm = T) 
  
  # For recrudescence, compare to previous episode only; returns NA if not estimate for previous
  C_unnorm <- sapply(2:n_epi, function(r) z[r-1, r, c("C")])
     
  # For relapse, compare to all preceding episodes... 
  L_unnorm <- apply(z[,-1,"L", drop = F], 2, max, na.rm = T) 
  
  # But if recrudescence has a high probability, overwrite all preceding with previous
  C_high <- C_unnorm > 0.5 & I_unnorm < 0.001
  L_prev <- sapply(2:n_epi, function(r) z[r-1, r, c("L")])
  L_unnorm[C_high] <- L_prev[C_high]
  
  # Compute normalised marginal posteriors; returns NA if recrudescence NA
  marg <- t(apply(cbind(C = C_unnorm, L = L_unnorm, I = I_unnorm), 1, function(x) x/sum(x)))
  
  # Check summation to one if not NA
  norm_check <- (rowSums(marg) - 1) <= .Machine$double.eps^0.5
  if (!all(norm_check, na.rm = T)) stop("Posterior probabilities do not sum to one")
  
  return(marg)
}


