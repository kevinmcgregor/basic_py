# Testing PY process function

source("~/research/pitman_yor/basic_py/py_functions.R")

# Number of samples
n.s <- 10000

# Concentration and discount values to try
conc <- c(0.01, 0.5, 1, 10)
dsct <- 0.2

# Sampling DP for fixed discount under different concentrations
samp <- vector("list", length=length(conc))
s1 <- matrix(0, n.s, length(conc))
for (i in 1:length(conc)) {
  cat("Concentration =", conc[i], "\n")
  samp[[i]] <- sample_py(n.s, conc[i], dsct)
}

par(mfrow=c(2,2))
for (i in 1:length(conc)) {
  plot(samp[[i]]$samp, pch=19, cex=0.5, col="blue",
       main=paste0("Concentration = ", conc[i], ", disc = ", dsct),
       xlab="Sample number", ylab="Species number")
}

# Number of customers at each table
par(mfrow=c(2,2))
for (i in 1:length(conc)) {
  plot(sort(samp[[i]]$s.c, decreasing = TRUE), type="l", lwd=3, col="blue",
       main=paste0("Concentration = ", conc[i], ", disc = ", dsct),
       xlab="Table (ordered by freq)", ylab="Number at table")
}

# Testing PY with specified frequency dist
fr <- c(0.1, 0.2, 0.05, 0.05, 0.3, 0.2, 0.1)
test.fr <- sample_py(1000, 0.5, 0.5, fr)
test.fr$s.c/sum(test.fr$s.c)

# Testing HPY sampler
n.top <- 1000
m <- 10
n <- rep(n.top, m)
conc.top <- 0.5
dsct.top <- 0.5
# Looping over different levels for lower-level params
conc.vals <- c(0.1, 1, 10, 100)
dsct.vals <- c(0.1, 0.5, 0.9)

lc <- length(conc.vals)
ld <- length(dsct.vals)

test.hpy <- vector("list", length=lc)
test.hpy <- lapply(test.hpy, function(x){return(vector("list", length=ld))})
pop.props <- vector("list", length=lc)
pop.props <- lapply(pop.props, function(x){return(vector("list", length=ld))})

for (i in 1:lc) {
  for (j in 1:ld) {
    cat("i =", i, ", j =", j, "\n")
    conc <- rep(conc.vals[i], m)
    dsct <- rep(dsct.vals[j], m)
    
    test.hpy[[i]][[j]] <- sample_hpy(m, n.top, n, conc.top, dsct.top, conc, dsct)
    
    # Comparing species distributions among top-level sample, and populations
    top.prop <- test.hpy[[i]][[j]]$top$s.c/sum(test.hpy[[i]][[j]]$top$s.c)
    pop.props[[i]][[j]] <- matrix(0, length(top.prop), m)
    for (k in 1:m) {
      pop.props[[i]][[j]][,k] <- test.hpy[[i]][[j]]$pop[[k]]$s.c/sum(test.hpy[[i]][[j]]$pop[[k]]$s.c)
    }
    pop.props[[i]][[j]] <- cbind(top.prop, pop.props[[i]][[j]])
  }
}

library(RColorBrewer)
par(mfrow=c(2,2))
d.ind <- 1

col <- brewer.pal(9, "Set1")
for (i in 1:lc) {
  barplot(pop.props[[i]][[d.ind]], col = col, names.arg = c("Top", paste0("Pop ", 1:m)),
          main = paste0("Conc = ", conc.vals[i], ", Dsct = ", dsct.vals[d.ind]))
}

