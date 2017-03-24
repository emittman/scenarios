s1 <- readRDS("small_circle_samples_sb.rds")
s2 <- readRDS("small_circle_samples_sd.rds")
t <- readRDS("small_circle_truth.rds")

source("constants.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(coda)

t$mles <- with(t, data$xty/(data$N/data$V))
t$est_prec <- with(t, sapply(1:G, function(g){
  data$N / (data$yty[g] + sum(mles[,g]^2) * data$N/data$V - 2 * t(mles[,g]) %*% data$xty[,g])
}))

##Traceplots of model parameters

#stick-breaking
mcmcs <- lapply(1:length(idx_save), function(g) mcmc(s1[[1]]$beta[1,,g]))
truths <- with(t, truth[1,zeta+1])
opar = par(mfrow=c(3,4))
for(g in 1:12){
  traceplot(mcmcs[[g]])
  abline(h=t$mles[1, idx_save[g]+1], col="red")
  abline(h=truths[g], col="blue", lty="dashed")
}
gpar <- par(opar) # reset plot grid
opar <- par(gpar)

mcmcs <- lapply(1:length(idx_save), function(g) mcmc(s1[[1]]$tau2[g,]))
for(g in 1:12){
  traceplot(mcmcs[[g]])
  abline(h=t$est_prec[idx_save[g]+1], col="red")
  abline(h=1/sigma_e^2, col="blue", lty=2)
}

gpar <- par(opar)
opar <- par(gpar)
#symmetric Dirichlet
mcmcs <- lapply(1:length(idx_save), function(g) mcmc(s2[[1]]$beta[1,,g]))
truths <- with(t, truth[1,zeta+1])
opar = par(mfrow=c(3,4))
for(g in 1:12){
  traceplot(mcmcs[[g]])
  abline(h=t$mles[1, idx_save[g]+1], col="red")
  abline(h=truths[g], col="blue", lty="dashed")
}
gpar <- par(opar) # reset plot grid
opar <- par(gpar)

mcmcs <- lapply(1:length(idx_save), function(g) mcmc(s2[[1]]$tau2[g,]))
for(g in 1:12){
  traceplot(mcmcs[[g]])
  abline(h=t$est_prec[idx_save[g]+1], col="red")
  abline(h=1/sigma_e^2, col="blue", lty=2)
}

gpar <- par(opar)


g_ran <- sample(length(idx_save), 1)
opar <- par(mfrow=c(2,1))
plot(t(s[[1]]$beta[,1:5000,g_ran]))
with(t, points(t(truth[,zeta[g_ran]+1]), col="blue"))
points(t(t$mles[,idx_save[g_ran]+1]), col="red")
plot(t(s[[1]]$beta[,5001:10000,g_ran]))
with(t, points(t(truth[,zeta[g_ran]+1]), col="blue"))
points(t(t$mles[,idx_save[g_ran]+1]), col="red")
ppar <- par(opar)

apply(s[[1]]$beta[,5001:10000,g_ran], 1, function(b) quantile(b, c(.01,.5,.99)))
with(t, truth[,zeta[g_ran]+1])
quantile(s[[1]]$tau2[g_ran, 1:2000], c(.01, .5, .99))
data.frame(beta=s[[1]]$beta[1,5001:10000,g_ran], tau2=s[[1]]$tau2[g_ran,5001:10000]) %>%
  ggplot(aes(beta, tau2)) + geom_hex(bins=30) +
  geom_point(data=data.frame(b0=0, t0=1), inherit.aes = F, aes(x=b0,y=t0), color="blue")+
  geom_point(data=data.frame(b0=t$mles[1,g_ran], t0=t$est_prec[g_ran]), inherit.aes=F,
             aes(x=b0,y=t0), color="red")

# weighted histograms
P <- with(s[[1]], data.frame(x1=as.numeric(P[,2,251:500]),
                             x2=as.numeric(P[,3,251:500]),
                             w=as.numeric(P[,1,251:500])))
true_means <- data.frame(x0=0, y0=0, r=5)

P %>%
ggplot(aes(x1, x2, weight=w)) + geom_hex(bins=25) + scale_fill_continuous()


#occupancy
opar <- par(ppar)
plot(prop.table(table(s[[1]]$max_id[5001:10000]+1)), main="max id", ylab="")
plot(prop.table(table(s[[1]]$num_occupied[5001:1000])), main="occupied clusters", ylab="")
ppar <- par(opar)
#check_probs
probs.df <- plyr::ldply(1:25, function(g){
  extern = mean(s[[1]]$beta[1,,g]>0 & s[[1]]$beta[2,,g]>0)
  intern = s[[2]]$probs[g]
  data.frame(extern=extern, intern=intern)
})

with(probs.df, all.equal(extern, intern))

#means
mean.mle <- data.frame(g = 1:25, group="mle", t(t$mles))
names(mean.mle)[3:4] <- c("beta1", "beta2")
mean.mle <- gather(mean.mle, key="component", value="value", 3:4)

mean.mcmc <- data.frame(g=1:25, group="mcmc", t(s[[2]]$means))
names(mean.mcmc)[3:4] <- c("beta1", "beta2")
mean.mcmc <- gather(mean.mcmc, key="component", value="value", 3:4)

mean.truth <- data.frame(g=1:25, t(t$truth[,t$zeta+1]))
names(mean.truth)[2:3] <- c("beta1", "beta2")
mean.truth <- gather(mean.truth, key="component", value="true_value", 2:3)

mean.df <- merge(mean.truth, rbind(mean.mle, mean.mcmc), by=c("g", "component"))

ggplot(mean.df, aes(x=factor(round(true_value, 1)), y=value, fill=group)) + 
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(breaks=c(-5,0,5)) +
  facet_wrap(~component)
  ggplot(mean.df, aes(x=beta1.truth, y=beta1.mcmc)) + geom_point()

#normal approximation
approx.mcmc <- with(s[[2]], plyr::ldply(1:25, function(g){
  mean <- means[,g]
  vars <- meansquares[,g] - means[,g]^2
  df <- data.frame(g=g, component=c("beta1","beta2"), mean = mean, sd = sqrt(vars))
}))

pdf("beta1_posteriors.pdf")
par(mfrow=c(5,5),
    oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(1, 1, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0))    # axis label at 2 rows distance, tick labels at 1 row

for(i in 1:25){
  hist(s[[1]]$beta[1,5001:10000,i], 100, prob=T, axes=F, main="")
  with(filter(approx.mcmc, g==i & component=="beta1"),
       curve(dnorm(x, mean, sd), add=T, lty=2, col=3))
  axis(side=1, at=-6:6, cex.axis=.7)
}
title(xlab = "value",
      ylab = "density",
      outer = TRUE, line = 3)
dev.off()

pdf("beta2_posteriors.pdf")
par(mfrow=c(5,5),
    oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(1, 1, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0))    # axis label at 2 rows distance, tick labels at 1 row

for(i in 1:25){
  hist(s[[1]]$beta[2,5001:10000,i], 100, prob=T, axes=F, main="")
  with(filter(approx.mcmc, g==i & component=="beta2"),
       curve(dnorm(x, mean, sd), add=T, lty=2, col=3))
  axis(side=1, at=-6:6, cex.axis=.7)
}
title(xlab = "value",
      ylab = "density",
      outer = TRUE, line = 3)
dev.off()