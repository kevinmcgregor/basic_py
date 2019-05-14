# Testing Dirichlet process function

source("~/research/pitman_yor/basic_py/py_functions.R")

# Number of samples
n.s <- 1000

# alpha values to try
alpha <- c(0.01, 1, 5, 10, 50, 100)

# Sampling DP from different alpha values
s <- matrix(0, n.s, length(alpha))
for (i in 1:length(alpha)) {
  
  s[,i] <- sample_dp(1, n.s, alpha[i])$samp
}

# Plotting samples over different alpha values
par(mfrow=c(3,2), mar=c(1,1,1,1))
for (i in 1:length(alpha)) {
  hist(s[,i], breaks=length(s[,i]), xlim=range(s),
       main=bquote(alpha ~ "=" ~ .(alpha[i])),
       xlab="Values")
}



