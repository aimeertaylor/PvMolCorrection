match_counting <- function(x){
  n_epi <- length(x)
  if(n_epi > 1) {
    
    markers <- names(x[[1]])
    n_markers <- length(markers)
    
    if (!all(sapply(x, function(x_epi) all(names(x_epi) == markers)))) {
      stop("Not all episodes typed at the same markers")}
    
    z <- sapply(2:n_epi, function(i){
      
      z_all <- sapply(1:(i-1), function(j) { # Compare to all previous
        mat <- sapply(1:n_markers, function(m) any(x[[i]][[m]] %in% x[[j]][[m]])) # get matches
        mat[is.na(x[[i]]) | is.na(x[[j]])] <- NA # replace NA matches with NA
        p <- mean(mat, na.rm = T) # proportion matched
        n <- sum(!is.na(mat)) # number of markers on which p is based
        return(rbind(p, n))
      })
      
      if (ncol(z_all) == 1 | all(is.na(z_all[1,]))) { # Only one or zero
        z_max <- z_all[,1, drop = FALSE] # Return first
      } else if (length(unique(z_all[1,])) == 1 & length(unique(z_all[2,])) == 1) { # All the same p and n
        z_max <- z_all[,1, drop = FALSE] # Return first 
      } else if (length(unique(z_all[1,])) == 1 & length(unique(z_all[2,])) > 1) { # All the same p 
        z_max <- z_all[, which.max(z_all[2,]), drop = FALSE] # Take maximum of n
      } else { # All different
        z_max <- z_all[, which.max(z_all[1,]), drop = FALSE] # Take maximum of p
      }
    
    return(z_max)
  
    })
    
    rownames(z) <- c("p", "n")
    if (!is.null(names(x))) colnames(z) <- names(x)[-1] # Return episode names if available
    return(z)
} else { z <- NA }
return(z)  
}

x1 <- list(initial = list(m1 = NA, m2 = NA), 
           recur = list(m1 = NA, m2 = NA), 
           recur2 = list(m1 = NA, m2 = NA))

x2 <- list(initial = list(m1 = c(1,2), m2 = c(1,4)), 
           recur = list(m1 = 1, m2 = 1))

x3 <- list(initial = list(m1 = c(0,2), m2 = c(1,4)), 
           recur = list(m1 = 1, m2 = 1))

x4 <- list(initial = list(m1 = c(0,1), m2 = 0), 
           recur = list(m1 = 2, m2 = c(1,2)))

match_counting(x1)
match_counting(x2)
match_counting(x3)
match_counting(x4)


