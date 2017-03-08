library(cudarpackage)
library(dplyr)

y <- readRDS("filtered_counts.rds")

log_yp1 <- log(y+1)

X <- matrix(c(1, 1, 0, 0,
              1, 0, 1, 1,
              1, -1, 0, 0,
              1, 0, 1, -1), 4, 4, byrow=T)
groups <- rep(1:4, each=4)

data <- formatData(counts=y, X=X, groups = groups)

ols <- lapply(1:nrow(log_yp1), function(g){
  fit <- lm(log_yp1[g,] ~ 0 + X[groups,])
  out <- c(coef(fit), sigma(fit))
  names(out) <- c("parent", "parent_hd", "hybrid", "hybrid_hd", "sigma")
  out  
})

ols <- do.call("rbind", ols)

#set priors for beta_k: multiply by 1.5 sd of ols, set prior mean to median of ols
mu_0 <- sapply(1:4, function(v) median(ols[,v]))
sigma_0 <- sapply(1:4, function(v) 1.5 *sd(ols[,v]) )

#set priors for tau_k: method of moments based on ols estimates of sigma_g
pr_tau2_var <- var(1/(ols[,'sigma']^2+.01)) #increase variance estimates slightly for stability
pr_tau2_mean <- mean(1/(ols[,'sigma']^2+.01))
a <- pr_tau2_mean^2 / pr_tau2_var
b <- pr_tau2_mean / pr_tau2_var

priors <- formatPriors(K=2500, prior_mean = mu_0, prior_sd = sigma_0, alpha = 20, a = a, b = b)
init_chain <- function(priors, G){
  beta <- with(priors, matrix(rnorm(V*K, rep(priors$mu_0, K), rep(sqrt(1/lambda2), K)), V, K))
  tau2 <- with(priors, rgamma(K, a, b))
  pi <- with(priors, rep(1/K, K))
  zeta <- with(priors, as.integer(sample(K, G, replace=T) - 1))
  formatChain(beta, pi, tau2, zeta)
}


chain <- init_chain(priors, data$G)
chain$C <- matrix(c(0, -1, 1, 1, # heterosis
                    0, -1, 1, -1,
                    0, 1, 1, 1,
                    0, 1, 1, -1), 4, 4, byrow=T)
chain$P <- as.integer(4) #nrow(C). Redundant, since P = V = 4 by default

idx_save <- sample(data$G, 20) - 1

out <- mcmc(data, priors, chain, n_iter = 1e4, n_save_P = 1e2, idx_save = idx_save, thin = 1, verbose = 0)
saveRDS(out, "output_paschold.rds")