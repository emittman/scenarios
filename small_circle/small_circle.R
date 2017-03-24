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
out1 <- mcmc(data, priors, methodPi="stickBreaking", n_iter = n_iter,
            n_save_P=1, idx_save = idx_save, thin = 1, verbose =0, C=C)
out2 <- mcmc(data, priors, methodPi="symmDirichlet", n_iter = n_iter,
             n_save_P=n_save_P, idx_save = idx_save, thin = 1, verbose = 0, C=C)
out3 <- mcmc(data, priors, methodPi="stickBreaking", alpha_fixed = F, n_iter = n_iter,
             n_save_P=1, idx_save = idx_save, thin = 1, verbose =0, C=C)
out4 <- mcmc(data, priors, methodPi="symmDirichlet", alpha_fixed = F, n_iter = n_iter,
            n_save_P=1, idx_save = idx_save, thin = 1, verbose =0, C=C)
saveRDS(out1, file="small_circle_samples_sb.rds")
saveRDS(out2, file="small_circle_samples_sd.rds")
saveRDS(list(truth = betas, zeta = zeta, data = data, priors = priors), file="small_circle_truth.rds")

