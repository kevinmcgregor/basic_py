# Functions for basic Pitman-Yor sampler


#' Sample from Dirichlet process
#'
#' @param n Integer specifying the number of sites to sample from
#' @param n.s Vector of integers specifying the number of individuals to sample from each site
#' @param b.dist Function that generates a sample from the base distribution
#' @param ... Additional parameters for the base distribution
#' @param alpha The concentration parameter, must be a positive real number
#'
#' @return
#' @export
#'
#' @examples
sample_dp <- function(n, n.s, alpha, b.dist = rnorm, ...) {
  # Temp: sample from just one site to start with
  n <- 1
  
  # Vector to contained sampled values
  samp <- rep(0, n.s)
  tot.spec <- 0
  
  for (i in 1:n.s) {
    #if (i%%100==0) cat("i =", i, "\n")
    
    samp.new <- rbinom(1, 1, alpha/(i-1+alpha))
    
    if (samp.new==1) {
      # Sample new value from base dist
      samp[i] <- b.dist(1, ...)
      tot.spec <- tot.spec + 1
    } else {
      # Sample one of the previous values
      samp[i] <- sample(samp[1:(i-1)], 1)
    }
  }
  
  return(list(samp=samp, tot.spec=tot.spec))
}

#' Sample from Pitman-Yor process
#'
#'
#' @param n.s Number of samples
#' @param conc Concentration parameter (must be greater than minus discount parameter)
#' @param dsct Discount parameter (must be between 0 and 1)
#' @param freq Frequency distribution from which to sample new species. Must be a vector of probabilities summing to 1
#'
#' @return
#' @export
#'
#' @examples
sample_py <- function(n.s, conc, dsct, freq=NULL) {
  if (!n.s>0 | !floor(n.s)==n.s) stop("n.s must be a positive integer")
  if (dsct<0 | dsct>1) stop("dsct must be between 0 and 1")
  if (conc<=(-dsct)) stop("conc must be greater than -dsct")
  if (!is.null(freq) & ((any(freq<0 | freq>1) | sum(freq)!=1))) 
    stop("freq must be a vector of probabilities with sum = 1")
  
  # Vector to contained sampled values
  samp <- rep(0, n.s)
  #Container for species counts
  if (is.null(freq)) {
    s.c <- NULL
  } else {
    s.c <- rep(0, length(freq))
  }
  t.c <- NULL #Table counts
  t.t <- 0 #Total tables
  t.s <- 0 #Total species
  
  for (i in 1:n.s) {
    p.new <- (conc+dsct*t.t)/(i-1+conc)
    samp.new <- rbinom(1, 1, p.new)
    
    if (samp.new==1) {
      # Sample new table
      if (is.null(freq)) {
        samp[i] <- t.s + 1
        s.c <- c(s.c, 1)
        t.c <- rbind(t.c, c(t.t+1, 1))
        t.s <- t.s + 1
      } else {
        # Sample species from top-level dist for new table
        mn <- c(rmultinom(1, 1, freq))
        wh.mn <- which(mn==1)
        samp[i] <- wh.mn
        t.s <- ifelse(s.c[wh.mn]>0, t.s, t.s + 1)
        s.c <- s.c + mn
        t.c <- rbind(t.c, c(wh.mn, 1))
      }
      t.t <- t.t + 1
    } else {
      # Sample existing table
      p <- pmax(0, (t.c[,2]-dsct)/(i-1-dsct*t.t))
      mn <- c(rmultinom(1, 1, p))
      wh.mn <- which(mn==1)
      samp[i] <- t.c[wh.mn, 1]
      s.c[t.c[wh.mn,1]] <- s.c[t.c[wh.mn,1]] + 1
      t.c[wh.mn, 2] <- t.c[wh.mn, 2] + 1
    }
  }
  
  colnames(t.c) <- c("Species", "Count")
  return(list(samp=samp, s.c=s.c, t.c=t.c, t.s=t.s, t.t=t.t))
}

#' Sample from hierarchical Pitman-Yor process
#'
#' @param m Number of populations (positive integer)
#' @param n.top Number of top-level samples (positive integer)
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
sample_hpy <- function(m, n.top, n, conc.top, dsct.top, conc, dsct, quiet=FALSE) {
  if (!m>0 | !floor(m)==m) stop("m must be a positive integer")
  if (!n.top>0 | !floor(n.top)==n.top) stop("n.top must be a positive integer")
  if (!all(n>0) | !all(floor(n)==n) | length(n)!=m) stop("n must be a vector of length m containing positive integers")
  if (dsct.top<0 | dsct.top>1) stop("dsct.top must be between 0 and 1")
  if (conc.top<=(-dsct.top)) stop("conc.top must be greater than -dsct.top")
  if (any(dsct<0 | dsct>1)) stop("all elements of dsct must be between 0 and 1")
  if (any(conc<=(-dsct))) stop("each element of conc must be less than or equal to the corresponding element of -dsct")
  
  # Sample top-level PY
  if (!quiet) cat("Sampling top-level PY", "\n")
  top.samp <- sample_py(n.top, conc.top, dsct.top)
  
  # Sample populations
  pop.samp <- vector("list", length=m)
  for (i in 1:m) {
    if (!quiet) cat("Sampling population", i, "\n")
    pop.samp[[i]] <- sample_py(n[m], conc[m], dsct[m], top.samp$s.c/sum(top.samp$s.c))
  }

  return(list(top=top.samp, pop=pop.samp))
}

