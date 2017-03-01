library(cudarpackage)

set.seed(22317)

G <- 100

N <- 10

sigma0 <- 5

sigma_e <- 1

betas <- sapply(1:G, function(g) c(sigma0 * cos(g*2*pi/G), sigma0 * sin(g*2*pi/G)))

X <- kronecker(diag(2), rep(1, reps))

y <- t(apply(betas, 2, function(b) rnorm(2*reps, X%*%b, sigma_e)))

data <- formatData(y, X, transform_y=identity)
priors <- formatPriors(K = G, prior_mean = rep(0,2), prior_sd = rep(sigma0,2), alpha = 5, a = 3, b = 2)
idx_save = 0:9 * 10
out <- mcmc(data, priors, n_iter = 1000, idx_save = idx_save, thin = 1, verbose = 0)

saveRDS(out, file="circle_samples.rds")
saveRDS(list(truth = betas, data = data, priors = priors), file="circle_truth.rds")