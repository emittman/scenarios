library(cudarpackage)

set.seed(22317)

source("constants.R")

betas <- sapply(1:K, function(g) c(sigma0 * cos(g*2*pi/K), sigma0 * sin(g*2*pi/K)))

X <- kronecker(diag(2), rep(1, NperV))

zeta <- sample(0:(K-1), G, replace=T)

y <- t(apply(betas[,zeta+1], 2, function(b) rnorm(2*NperV, X%*%b, sigma_e)))

init_chain <- function(priors, G){
  beta <- with(priors, matrix(rnorm(V*K, rep(priors$mu_0, K), rep(sqrt(1/lambda2), K)), V, K))
  tau2 <- with(priors, rgamma(K, a, b))
  pi <- with(priors, rep(K, 1/K))
  zeta <- with(priors, as.integer(sample(K, G, replace=T) - 1))
  formatChain(beta, pi, tau2, zeta)
}

data <- formatData(y, X, transform_y=identity)
priors <- formatPriors(K = G, prior_mean = rep(0,2), prior_sd = rep(sigma0,2), alpha = 1, a = 3, b = 2)
chain <- init_chain(priors, G)
out <- mcmc(data, priors, chain, n_iter = 10000, n_save_P, idx_save = idx_save, thin = 1, verbose = 0)

saveRDS(out, file="small_circle_samples.rds")
saveRDS(list(truth = betas, data = data, priors = priors), file="small_circle_truth.rds")

