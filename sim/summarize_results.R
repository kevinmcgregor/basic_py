# Read in results from all simulations and summarize results
library(ggplot2)

in.dir <- "/mnt/GREENWOOD_SCRATCH/kevin.mcgregor/pitman_yor/basic_py/sim"

# Values of n and concentration parameter to loop over
# (input filenames are based on these)
n.vals <- c(10, 100, 1000)
conc.vals <- c(0.1, 1, 10, 50)
dsct.vals <- 0.5
n.rep <- 100 # Number of replications of each sim

n.n <- length(n.vals)
n.c <- length(conc.vals)

# Containers
bias.top <- bias.local <- matrix(0, n.rep, n.n*n.c)
bias.top.exact <- bias.local.exact <- matrix(0, n.rep, n.n*n.c)

for (nn in 1:n.n) {
  cat("n =", n.vals[nn], "\n")
  for (cc in 1:n.c) {
    cat("conc =", conc.vals[cc], "\n")
    fname <- paste0(in.dir, "/", "bias_n", n.vals[nn], "_conc", conc.vals[cc], "_dsct", 
                    dsct.vals, ".RData")
    load(fname)
    
    # Calculating bias
    ind <- (nn-1)*n.c+cc
    bias.top[,ind] <- (obs.n.sp - expect.n.sp)/expect.n.sp
    bias.local[,ind] <- (obs.n.sp.lower - expect.n.sp.lower)/expect.n.sp.lower
    bias.top.exact[,ind] <- (obs.n.sp - expect.n.sp.exact)/expect.n.sp.exact
    bias.local.exact[,ind] <- (obs.n.sp.lower - expect.n.sp.lower.exact)/expect.n.sp.lower.exact
  }
}

# Setting up data frames for ggplot
f.n <- factor(n.vals, levels=c("10", "100", "1000"))
f.c <- factor(conc.vals, levels=c("0.1", "1", "10", "50"))

local.df <- data.frame(Bias=c(bias.local),
                     n=rep(f.n, each=n.rep*n.c),
                     Concentration=rep(rep(f.c, each=n.rep), n.n))

local.df.exact <- local.df
local.df.exact$Bias <- c(bias.local.exact)



plot.title.local = paste0("Local-level bias: approximate EV, ", "Discount = ", dsct.vals[1])
local.plot = ggplot(local.df, aes(n, Bias, fill=Concentration)) + geom_boxplot() + #facet_wrap( ~ z*n.otu, ncol=2) + 
  theme_bw() +
  theme(legend.position = "bottom",text = element_text(size=20), axis.text.x = element_text(size=16,angle=0,vjust=1,hjust=0.5)) + 
  labs(title=plot.title.local) +
  xlab("Number local samples") + ylab("Bias") + cowplot::panel_border()

plot.title.local.exact = paste0("Local-level bias: exact EV", "Discount = ", dsct.vals[1])
local.plot.exact = ggplot(local.df.exact, aes(n, Bias, fill=Concentration)) + geom_boxplot() + #facet_wrap( ~ z*n.otu, ncol=2) + 
  theme_bw() +
  theme(legend.position = "bottom",text = element_text(size=20), axis.text.x = element_text(size=16,angle=0,vjust=1,hjust=0.5)) + 
  labs(title=plot.title.local.exact) +
  xlab("Number local samples") + ylab("Bias") + cowplot::panel_border()


