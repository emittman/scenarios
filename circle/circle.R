library(cudarpackage)

set.seed(22317)

G <- 500

reps <- 4

sigma0 <- 5

sigma_e <- 1

betas <- sapply(1:G, function(g) c(sigma0 * cos(g*2*pi/G), sigma0 * sin(g*2*pi/G)))

X <- kronecker(diag(2), rep(1, reps))

y <- t(apply(betas, 2, function(b) rnorm(2*reps, X%*%b, sigma_e)))

data <- formatData(y, X, transform_y=identity)
priors <- formatPriors(K = 200, prior_mean = rep(0,2), prior_sd = rep(sigma0,2), alpha = 5, a = 2, b = 2)

out <- mcmc(data, priors, n_iter = 10000, idx_save = 0:(G-1), thin = 10, verbose = 0)

saveRDS(out, file="circle_samples.rds")
saveRDS(list(truth = betas, data = data, priors = priors), file="circle_truth.rds")