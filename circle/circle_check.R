s <- readRDS("circle_samples.rds")
t <- readRDS("circle_truth.rds")
library(coda)
mcmcs <- lapply(1:10*20, function(g) mcmc(s$beta[1,,g]))
truths <- t$truth[1, 1:10*20]
par(mfrow=c(3,4))
for(g in 1:10){
  traceplot(mcmcs[[g]], ylim=c(-15,15))
  abline(h=truths[g], col="red")
}
dev.off()

g_ran <- sample(200, 1)
par(mfrow=c(2,1))
plot(t(s$beta[,1:500,g_ran]))
points(t(t$truth[,g_ran]), col="red")
plot(t(s$beta[,501:1000,g_ran]))
points(t(t$truth[,g_ran]), col="red")
dev.off()

sd(s$beta[,501:1000,g_ran])

plot(s$beta[1,501:1000,g_ran], s$tau2[g_ran,501:1000])
