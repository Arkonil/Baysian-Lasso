rm(list=ls())
set.seed(0)

library(MASS)
source("./src/original_lasso.R")
source("./src/fused_lasso.R")
source("./src/elasticnet_lasso.R")

compute.R1 = function(p) {
    R1 = matrix(, nrow = p, ncol = p)
    for (i in 1:p) {
        for (j in 1:p) {
            R1[i, j] = 0.5 ^ abs(i-j)
        }
    }
    return(R1)
}

true.beta = c(3, 1.5, 0, 0, 2, 0, 0, 0)
true.sigmaSq = 9

p = length(true.beta)
R = compute.R1(p)

n = 200
X = mvrnorm(n, mu=rep(0, p), Sigma=R)

N = 50

result = list(
    original = gibbs.original.lasso(X, true.beta, true.sigmaSq, N),
    elasticnet = gibbs.elasticnet.lasso(X, true.beta, true.sigmaSq, N),
    fused = gibbs.fused.lasso(X, true.beta, true.sigmaSq, N)
)

saveRDS(result, "example_1_results.RData")