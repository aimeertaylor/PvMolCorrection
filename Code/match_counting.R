match_counting <- function(x){
  n_epi <- length(x)
  if(n_epi > 1) {
    
    markers <- names(x[[1]])
    n_markers <- length(markers)
    
    if (!all(sapply(x, function(x_epi) all(names(x_epi) == markers)))) {
      stop("Not all episodes typed at the same markers")}
    
    z <- sapply(1:(n_epi-1), function(i){
      m <- sapply(1:n_markers, function(j) any(x[[i+1]][[j]] %in% x[[i]][[j]])) # get matches
      m[is.na(x[[1]]) | is.na(x[[2]])] <- NA # replace NA matches with NA
      p <- mean(m, na.rm = T) # proportion of matches
      n <- sum(!is.na(m)) # number of markers with comparable data
      return(rbind(p, n))
    })
    rownames(z) <- c("p", "n")
    if (!is.null(names(x))) colnames(z) <- names(x)[-1] # Return episode names if available
  } else { z <- NA }
  return(z)  
}

x1 <- list(initial = list(m1 = NA, m2 = NA), 
          recur = list(m1 = NA, m2 = NA))

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


