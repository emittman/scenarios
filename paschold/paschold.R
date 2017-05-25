library(cudarpackage)
library(dplyr)
set.seed(1001001001)
y <- readRDS("filtered_counts.rds")

log_yp1 <- log(y+1)

X <- matrix(c(1, 1, 0, 0,
              1, 0, 1, 1,
              1, -1, 0, 0,
              1, 0, 1, -1), 4, 4, byrow=T)
groups <- rep(1:4, each=4)

data <- formatData(counts=y, X=X, groups = groups)

estimates <- indEstimates(data)

priors <- formatPriors(K=5000, estimates = estimates, A=1.5, B=1.5/sqrt(nrow(y)))

C <- list(hph = matrix(c(0, -1, 1, 0, 
                         0, 1,  1, 0),2,4, byrow=T),
          lph = matrix(c(0, -1, -1, 0,
                         0, 1, -1, 0),2,4, byrow=T),
          hybpos = matrix(c(0, 0, 1, 0),1,4, byrow=T),
          hybneg = matrix(c(0, 0, -1, 0),1, 4, byrow=T)) 

idx_save <- 1:10*200-1
n_iter <- 5000
warmup <- 1000
n_save_P <- 300
s <- mcmc(data = data, priors = priors, methodPi="stickBreaking", n_iter = n_iter, C = C,
                        alpha_fixed = FALSE, n_save_P = n_save_P, warmup= warmup, idx_save = idx_save,
                        thin = 10, verbose = 0, estimates=estimates)

saveRDS(s, "samplesSB_524.rds")
