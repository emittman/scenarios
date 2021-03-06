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
priors <- formatPriors(K = G, prior_mean = rep(0,2), prior_sd = rep(sigma0,2), a = 3, b = 2, A = 1, B = 1/sqrt(G))
control <- formatControl(n_iter=n_iter, thin=10, warmup=warmup, methodPi="stickBreaking",
                         idx_save=idx_save, n_save_P=n_save_P, alpha_fixed = TRUE)

saveRDS(list(truth = betas, zeta = zeta, data = data, priors = priors), file="small_circle_truth.rds")

out1 <- mcmc(data, priors, control, verbose =0, C=C, alpha=10)
control$methodPi="symmDirichlet"
out2 <- mcmc(data, priors, control, verbose = 0, C=C, alpha=10)
saveRDS("out1", file="small_circle_samples_sb_fixed.rds")
saveRDS("out2", file="small_circle_samples_sd_fixed.rds")

control$alpha_fixed=FALSE
out4 <- mcmc(data, priors, control, verbose =0, C=C)
control$methodPi="stickBreaking"
out3 <- mcmc(data, priors, control, verbose =0, C=C)

saveRDS(out3, file="small_circle_samples_sb.rds")
saveRDS(out4, file="small_circle_samples_sd.rds")

