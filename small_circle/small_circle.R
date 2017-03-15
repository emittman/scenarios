library(cudarpackage)

set.seed(22317)

source("constants.R")

betas <- sapply(1:K, function(g) c(sigma0 * cos(g*2*pi/K), sigma0 * sin(g*2*pi/K)))

X <- kronecker(diag(2), rep(1, NperV))

zeta <- sample(0:(K-1), G, replace=T)

y <- t(apply(betas[,zeta+1], 2, function(b) rnorm(2*NperV, X%*%b, sigma_e)))

C <- list(matrix(c(1, 0), 1, 2),
          matrix(c(0, 1), 1, 2))

data <- formatData(y, X, transform_y=identity)
priors <- formatPriors(K = G, prior_mean = rep(0,2), prior_sd = rep(sigma0,2), alpha = 1, a = 3, b = 2)
out <- mcmc(data, priors, n_iter = 10000, n_save_P=n_save_P, idx_save = idx_save, thin = 1, verbose = 0, C=C)

saveRDS(out, file="small_circle_samples.rds")
saveRDS(list(truth = betas, zeta = zeta, data = data, priors = priors), file="small_circle_truth.rds")

