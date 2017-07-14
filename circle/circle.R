library(cudarpackage)

set.seed(22317)

source("constants.R")

betas <- sapply(1:G, function(g) c(sigma0 * cos(g*2*pi/G), sigma0 * sin(g*2*pi/G)))

X <- kronecker(diag(2), rep(1, NperV))

y <- t(apply(betas, 2, function(b) rnorm(2*NperV, X%*%b, sigma_e)))

data <- formatData(y, X, transform_y=identity)
priorsfixed <- formatPriors(K = G/2, prior_mean = rep(0,2), prior_sd = rep(sigma0,2),
                       a = 3, b = 2, A = 10000, B = 100)

system.time(out1 <- mcmc(data, priorsfixed, methodPi="stickBreaking", n_iter = n_iter, n_save_P = n_save_P,
            idx_save = idx_save, thin = thin, verbose = 0, warmup=warmup))

system.time(out2 <- mcmc(data, priorsfixed, methodPi="symmDirichlet", n_iter = n_iter, n_save_P = n_save_P,
                         idx_save = idx_save, thin = thin, verbose = 0, alpha = 10, warmup=warmup))

priors <- formatPriors(K = G/2, prior_mean = rep(0,2), prior_sd = rep(sigma0,2),
                            a = 3, b = 2, A = 10, B = .1)

system.time(out3 <- mcmc(data, priors, methodPi="stickBreaking", alpha_fixed=FALSE, n_iter = n_iter, n_save_P = n_save_P,
             idx_save = idx_save, thin = thin, verbose = 0, warmup=warmup))

system.time(out4 <- mcmc(data, priors, methodPi="symmDirichlet", alpha_fixed=FALSE, n_iter = n_iter, n_save_P = n_save_P,
             idx_save = idx_save, thin = thin, verbose = 0, warmup=warmup))


saveRDS(list(truth = betas, data = data, priors = priors), file="circle_truth.rds")
saveRDS(out1, file="circle_samples.rds")
saveRDS(out2, file="circle_samples_SD.rds")
saveRDS(out3, file="circle_samples_freealpha.rds")
saveRDS(out4, file="circle_samples_SD_freealpha.rds")
