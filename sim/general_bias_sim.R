# Code for general simulation for HPY sampler
library(parallel)
source("~/research/pitman_yor/basic_py/R/hpy_functions.R")

cmd.line <- TRUE
n.cores <- 10

# Reading command-line arguments
# (Filename with list of params)
if (cmd.line) {
  arg.file <- commandArgs(trailingOnly = TRUE)
} else {
  arg.file <- "/mnt/GREENWOOD_BACKUP/home/kevin.mcgregor/research/pitman_yor/basic_py/sim/arg_files/af1.txt"
}

args <- read.table(arg.file, header = FALSE, stringsAsFactors = FALSE)$V1
n.sim <- as.numeric(args[1])
m <- as.numeric(args[2])
n <- as.numeric(args[3])
conc <- as.numeric(args[4])
dsct <- as.numeric(args[5])
out.dir <- args[6]

# Run HPY simulations in parallel
sim.results <- mclapply(1:n.sim, sample_hpy, m=m, n=rep(n, m),
                        conc.top=conc, dsct.top=dsct,
                        conc=rep(conc,m), dsct=rep(dsct,m),
                        mc.cores = n.cores)

# Observed numbers of species at top level
obs.n.sp <- unlist(lapply(sim.results, function(x){x$n.sp}))

# Theoretical means for number of tables at top level (using approximate formula)
expect.n.sp <- unlist(lapply(sim.results, function(x) {expected_tables(x$n.samp, conc, dsct)}))
# Theoretical means for number of tables at top level (using exact formula)
expect.n.sp.exact <- unlist(lapply(sim.results, function(x) {
  expected_tables(x$n.samp, conc, dsct, method="exact")
}))

# Same thing for local level
obs.n.sp.lower <- unlist(lapply(sim.results, function(x){mean(x$t.tab)}))
expect.n.sp.lower <- expected_tables(n, conc, dsct)
expect.n.sp.lower.exact <- expected_tables(n, conc, dsct, method="exact")

save(m, n, conc, dsct, obs.n.sp, obs.n.sp.lower, expect.n.sp, expect.n.sp.exact,
     expect.n.sp.lower, expect.n.sp.lower.exact,
     file=paste0(out.dir, "/", "bias_n", n, "_conc", conc, "_dsct", dsct, ".RData"))
