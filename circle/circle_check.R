s <- readRDS("circle_samples_freealpha.rds")
t <- readRDS("circle_truth.rds")

source("constants.R")

library(coda)

t$mles <- with(t, data$xty/(data$N/data$V))
t$est_prec <- with(t, sapply(1:G, function(g){
  data$N / (data$yty[g] + sum(mles[,g]^2) * data$N/data$V - 2 * t(mles[,g]) %*% data$xty[,g])
}))
mcmcs <- lapply(1:length(idx_save), function(g) mcmc(s[[1]]$beta[1,,g]))
truths <- t$truth[1, idx_save+1]
opar = par(mfrow=c(3,4))
for(g in 1:10){
  traceplot(mcmcs[[g]])
  abline(h=t$mles[1, idx_save[g]+1], col="red")
  abline(h=truths[g], col="blue", lty="dashed")
}
gpar <- par(opar)

mcmcs <- lapply(1:length(idx_save), function(g) mcmc(s[[1]]$tau2[g,]))
opar = par(gpar)
for(g in 1:10){
  traceplot(mcmcs[[g]])
  abline(h=t$est_prec[idx_save[g]+1], col="red")
  abline(h=1/sigma_e^2, col="blue", lty=2)
}
par(opar)

g_ran <- sample(length(idx_save), 1)
par(mfrow=c(2,1))
plot(t(s[[1]]$beta[,1:1000,g_ran]))
points(t(t$truth[,idx_save[g_ran]+1]), col="blue")
points(t(t$mles[,idx_save[g_ran]+1]), col="red")
plot(t(s[[1]]$beta[,1000:2000,g_ran]))
points(t(t$truth[,idx_save[g_ran]+1]), col="blue")
points(t(t$mles[,idx_save[g_ran]+1]), col="red")
par(opar)

apply(s[[1]]$beta[,1001:2000,g_ran], 1, function(b) quantile(b, c(.01,.5,.99)))
t$truth[,idx_save[g_ran]+1]
quantile(s[[1]]$tau2[g_ran, 1:2000], c(.01, .5, .99))
plot(s[[1]]$beta[1,1001:2000,g_ran], s[[1]]$tau2[g_ran,1001:2000])

# weighted histograms
library(ggplot2)
library(dplyr)
library(ggforce)
df <- data.frame(x=runif(1e6))
df <- dplyr::mutate(df, y=runif(1e6, 0, x))
df <- dplyr::mutate(df, weight = dnorm(y, .3, .1))
ggplot(df, aes(x, y, weight=weight)) + geom_hex(bins=45) + 
  theme_bw() + coord_fixed()

beta_df <- data.frame(x = as.numeric(s[[1]]$beta[1,1001:2000,]),
                      y = as.numeric(s[[1]]$beta[2,1001:2000,]))
ggplot(beta_df, aes(x, y)) + geom_hex() + coord_fixed()

P <- with(s[[1]], data.frame(x1=as.numeric(P[,2,]),x2=as.numeric(P[,3,]), w=as.numeric(P[,1,])))
true_means <- data.frame(x0=0, y0=0, r=5)

P %>%
ggplot(aes(x1, x2, weight=w)) + geom_hex(bins=25) + scale_fill_continuous(trans="log", low="white", high="black") +
  geom_circle(data=true_means, inherit.aes = F, aes(x0=x0,y0=y0,r=r), color="red", size=1) + theme_bw()


#occupancy
hist(s[[1]]$max_id[1001:2000], 100)
hist(s[[1]]$num_occupied[1001:2000], 100)
