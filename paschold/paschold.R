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

estimates <- indEstimates(data)

priors <- formatPriors(K=5000, estimates = estimates, A=3, B=3/sqrt(nrow(y)))

C <- list(hph = matrix(c(0, -1, 1, 0, 
                         0, 1,  1, 0),2,4, byrow=T),
          lph = matrix(c(0, -1, -1, 0,
                         0, 1, -1, 0),2,4, byrow=T),
          hybpos = matrix(c(0, 0, 1, 0),1,2, byrow=T),
          hybneg = matrix(c(0, 0, -1, 0),1, 2, byrow=T)) 

idx_save <- 1:10*200-1
n_iter <- 30000
warmup <- 5000
n_save_P <- 300
out_symm <- list()
s <- mcmc(data = data, priors = priors, methodPi="stickBreaking", n_iter = n_iter, C = C,
                        alpha_fixed = FALSE, n_save_P = n_save_P, warmup= warmup, idx_save = idx_save,
                        thin = 10, verbose = 0, estimates=estimates)

saveRDS(out_symm, "samplesSB_522.rds")
