# Functions for hierarchical Pitman-Yor sampler
# This version treats the top-level as latent, so the number of
# samples in the top-level is a random variable.


#' Sample from hierarchical Pitman-Yor process
#'
#' @param m Number of populations (positive integer)
#' @param n Vector of length m containing number of samples in each population
#' @param conc.top Top-level concentration parameter (must be greater than -dsct.top)
#' @param dsct.top Top-level discount parameter (must be between 0 and 1)
#' @param conc Vector of length m containing concentration parameters for each population 
#' @param dsct Vector of length m containing discount parameters for each population

#'
#' @return
#' @export
#'
#' @examples
sample_hpy <- function(m, n, conc.top, dsct.top, conc, dsct, quiet=FALSE) {
  if (!m>0 | !floor(m)==m) stop("m must be a positive integer")
  if (!all(n>0) | !all(floor(n)==n) | length(n)!=m) stop("n must be a vector of length m containing positive integers")
  if (dsct.top<0 | dsct.top>1) stop("dsct.top must be between 0 and 1")
  if (conc.top<=(-dsct.top)) stop("conc.top must be greater than -dsct.top")
  if (any(dsct<0 | dsct>1)) stop("all elements of dsct must be between 0 and 1")
  if (any(conc<=(-dsct))) stop("each element of conc must be less than or equal to the corresponding element of -dsct")
  
  # Initializing containers
  n.s <- 0 # Total number of species in top-level process
  t.t <- rep(0, m) # Total number of tables in each population
  tab <- vector("list", length=m) # Container for tables in each population
  samp.top <- NULL # Samples for top-level process
  samp <- vector("list", length=m) # List to hold samples for each population
  
  # Loop through populations and samples in each population
  for (i in 1:m) {
    samp[[i]] <- rep(0, n[i])
    for (j in 1:n[i]) {
      p.new <- (conc[i]+dsct[i]*t.t[i])/(j-1+conc[i])
      samp.new <- rbinom(1, 1, p.new)
      
      if (samp.new==1) {
        s <- sample_top(samp.top, conc.top, dsct.top)
        tab[[i]] <- rbind(tab[[i]], c(s, 1))
        samp[[i]][j] <- s
        t.t[i] <- t.t[i] + 1
        if (s > n.s) {
          # New species at top level
          samp.top <- c(samp.top, 1)
          n.s <- n.s + 1
        } else {
          samp.top[s] <- samp.top[s] + 1
        }
      } else {
        ind <- sample_local(tab[[i]], dsct[i])
        samp[[i]][j] <- tab[[i]][ind,1]
        tab[[i]][ind,2] <- tab[[i]][ind,2] + 1
      }
    }
  }
  
  return(list(samp.pop=samp, samp.top=samp.top, tab=tab, t.t=t.t, n.s=n.s))
}

#' Sample from top-level process
#'
#' @param samp Vector of current top-level counts
#' @param conc Top-level concentration parameter
#' @param dsct Top-level discount parameter
#'
#' @return
#' @export
#'
#' @examples
sample_top <- function(samp, conc, dsct) {
  n.s <- sum(samp)
  t.t <- length(samp)
  p.new <- (conc+dsct*t.t)/(n.s+conc)
  samp.new <- rbinom(1, 1, p.new)
  
  #if(is.na(samp.new)) browser()
  
  if (samp.new==1) {
    r <- t.t + 1
  } else {
    r <- sample(t.t, 1, prob=samp/n.s)
  }
  
  return(r)
}

#' Sample local ancestor
#'
#' @param tab Table for this population
#' @param dsct Discount parameter for this population
#'
#' @return
#' @export
#'
#' @examples
sample_local <- function(tab, dsct) {
  n.s <- sum(tab[,2])
  t.t <- NROW(tab)
  p <- (tab[,2]-dsct)/(n.s-dsct*t.t)
  return(sample(t.t, 1, prob=p))
}

