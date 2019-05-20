# Small simulation for HPY to check if there's bias in # of tables
library(parallel)

source("~/research/pitman_yor/basic_py/R/hpy_functions.R")

n.sim <- 100

m <- 10
n <- 1000
conc <- 0.5
dsct <- 0.5

sim.results <- mclapply(1:n.sim, sample_hpy, m=m, n=rep(n, m),
                        conc.top=conc, dsct.top=dsct,
                        conc=rep(conc,m), dsct=rep(dsct,m),
                        mc.cores = 10)

# Observed numbers of species at top level
obs.n.sp <- unlist(lapply(sim.results, function(x){x$n.sp}))

# Theoretical means for number of tables at top level (using approximate formula)
expect.n.sp <- unlist(lapply(sim.results, function(x) {expected_tables(x$n.samp, conc, dsct)}))
# Theoretical means for number of tables at top level (using exact formula)
expect.n.sp.exact <- unlist(lapply(sim.results, function(x) {
  expected_tables(x$n.samp, conc, dsct, method="exact")
}))

# Plotting top-level results for approximate theoretical values
plot(expect.n.sp, obs.n.sp)
abline(0,1,col="red", lwd=3)

boxplot(obs.n.sp-expect.n.sp, col="lightblue")
abline(h=0, col="red", lwd=3)

# Plotting top-level results for exact theoretical values
plot(expect.n.sp.exact, obs.n.sp)
abline(0,1,col="red", lwd=3)

boxplot(obs.n.sp-expect.n.sp.exact, col="lightblue")
abline(h=0, col="red", lwd=3)


# Plotting local results for approximate theoretical values
obs.n.sp.lower <- unlist(lapply(sim.results, function(x){mean(x$t.tab)}))
expect.n.sp.lower <- expected_tables(n, conc, dsct)
expect.n.sp.lower.exact <- expected_tables(n, conc, dsct, method="exact")

plot(obs.n.sp.lower)
abline(h=expect.n.sp.lower,col="red", lwd=3)

boxplot(obs.n.sp.lower-expect.n.sp.lower, col="lightblue")
abline(h=0, col="red", lwd=3)

# Plotting top-level results for exact theoretical values
plot(obs.n.sp.lower)
abline(h=expect.n.sp.lower.exact,col="red", lwd=3)

boxplot(obs.n.sp.lower-expect.n.sp.lower.exact, col="lightblue")
abline(h=0, col="red", lwd=3)

# Comparing approximate vs. exact theoretical values for fixed N
n.fix <- 1000
dsct.fix <- 0.5
conc.vals <- seq(0.1, 20, by=0.1)
approx.ev <- rep(0, length(conc.vals))
exact.ev <- rep(0, length(conc.vals))
for (v in 1:length(conc.vals)) {
  approx.ev[v] <- expected_tables(n.fix, conc.vals[v], dsct.fix)
  exact.ev[v] <- expected_tables(n.fix, conc.vals[v], dsct.fix, method="exact")
}

rge.x <- range(exact.ev)
plot(exact.ev, approx.ev, ylim=rge.x, type="l", lwd=3)
abline(0, 1, lty=2, col="red", lwd=2)
axis(3, at=seq(rge.x[1], rge.x[2], length.out = 5),
     labels = conc.vals)
