s <- readRDS("circle_samples.rds")
t <- readRDS("circle_truth.rds")
library(coda)
mcmcs <- lapply(1:10*5, function(g) mcmc(s$beta[1,,g]))
truths <- t$truth[1, 1:10*5]
opar = par(mfrow=c(3,4))
for(g in 1:10){
  traceplot(mcmcs[[g]], ylim=c(-8,8))
  abline(h=truths[g], col="red")
}
par(opar)

mcmcs <- lapply(1:10*50, function(g) mcmc(s$tau2[g,]))
opar = par(mfrow=c(3,4))
for(g in 1:10){
  traceplot(mcmcs[[g]], ylim=c(0,2))
  abline(h=1, col="red")
}
par(opar)

g_ran <- sample(200, 1)
par(mfrow=c(2,1))
plot(t(s$beta[,1:500,g_ran]))
points(t(t$truth[,g_ran]), col="red")
plot(t(s$beta[,501:1000,g_ran]))
points(t(t$truth[,g_ran]), col="red")
par(opar)

apply(s$beta[,501:1000,g_ran], 1, function(b) quantile(b, c(.01,.5,.99)))
t$truth[,g_ran]
quantile(s$tau2[g_ran, 501:1000], c(.01, .5, .99))
plot(s$beta[1,501:1000,g_ran], 1/s$tau2[g_ran,501:1000])
