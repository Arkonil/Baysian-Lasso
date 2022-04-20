rm(list=ls())
set.seed(0)

library(MASS)
source("./src/original_lasso.R")
source("./src/fused_lasso.R")
source("./src/elasticnet_lasso.R")

compute.R2 = function(p) {
    R2 = matrix(0.5, nrow = p, ncol = p) + diag(0.5, p)
    return(R2)
}

true.beta = rep(c(0, 2), times=2, each=10)
true.sigmaSq = 15^2

p = length(true.beta)
R = compute.R2(p)

n = 500
X = mvrnorm(n, mu=rep(0, p), Sigma=R)

N = 50

result = list(
    original = gibbs.original.lasso(X, true.beta, true.sigmaSq, N),
    elasticnet = gibbs.elasticnet.lasso(X, true.beta, true.sigmaSq, N),
    fused = gibbs.fused.lasso(X, true.beta, true.sigmaSq, N)
)

saveRDS(result, "example_3_results.RData")