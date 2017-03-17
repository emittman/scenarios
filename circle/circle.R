library(cudarpackage)

set.seed(22317)

source("constants.R")

betas <- sapply(1:G, function(g) c(sigma0 * cos(g*2*pi/G), sigma0 * sin(g*2*pi/G)))

X <- kronecker(diag(2), rep(1, NperV))

y <- t(apply(betas, 2, function(b) rnorm(2*NperV, X%*%b, sigma_e)))

data <- formatData(y, X, transform_y=identity)
priors <- formatPriors(K = G, prior_mean = rep(0,2), prior_sd = rep(sigma0,2), alpha = 10, a = 3, b = 2)

out <- mcmc(data, priors, weightMethod="stickBreaking", n_iter = n_iter, n_save_P = n_save_P,
            idx_save = idx_save, thin = thin, verbose = 0)

out2 <- mcmc(data, priors, weightMethod="symmDirichlet", n_iter = n_iter, n_save_P = n_save_P,
             idx_save = idx_save, thin = thin, verbose = 2)

saveRDS(list(truth = betas, data = data, priors = priors), file="circle_truth.rds")
saveRDS(out, file="circle_samples.rds")
saveRDS(out2, file="circle_samples_SD.rds")