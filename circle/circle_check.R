s <- readRDS("circle_samples.rds")

library(coda)
mcmcs <- lapply(1:10, function(g) mcmc(s$beta[1,,g]))
par(mfrow=c(3,4))
for(g in 1:10){
  traceplot(mcmcs[[g]], ylim=c(-10,10))
}
dev.off()
