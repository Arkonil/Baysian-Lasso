library(MASS)
library(statmod)
library(IRdisplay)

calculate.Sigma.beta.inv = function(tauSq.inv, omegaSq.inv) {
    p = length(tauSq.inv)
    
    # Constructing Σβ inverse matrix
    Sigma.beta.inv = diag(tauSq.inv)
    Sigma.beta.inv[1:p-1, 1:p-1] = Sigma.beta.inv[1:p-1, 1:p-1] + diag(omegaSq.inv)
    Sigma.beta.inv[2:p, 2:p] = Sigma.beta.inv[2:p, 2:p] + diag(omegaSq.inv)
    Sigma.beta.inv[1:p-1, 2:p] = Sigma.beta.inv[1:p-1, 2:p] - diag(omegaSq.inv)
    Sigma.beta.inv[2:p, 1:p-1] = Sigma.beta.inv[2:p, 1:p-1] - diag(omegaSq.inv)
    return(Sigma.beta.inv)
}

fused.beta.update = function(sigmaSq, tauSq.inv, omegaSq.inv, X, y, sample.size=1) {
    p = length(tauSq.inv)
    
    # Constructing Σβ inverse matrix
    Sigma.beta.inv = calculate.Sigma.beta.inv(tauSq.inv, omegaSq.inv)
    
    A = solve(t(X) %*% X + Sigma.beta.inv)

    mu0 = A %*% t(X) %*% y
    sigmaSq0 = sigmaSq * A

    beta = mvrnorm(sample.size, mu=mu0, Sigma=sigmaSq0)
    return(beta)
}

fused.tauSq.inv.update = function(beta, sigmaSq, lambda1Sq) {
    p = length(beta)

    means = sqrt((lambda1Sq * sigmaSq) / (beta^2))
    shape = lambda1Sq

    tauSq.inv = numeric(p)
    for (i in 1:p) {
        tauSq.inv[i] = rinvgauss(1, means[i], shape)
    }

    return(tauSq.inv)
}

fused.omegaSq.inv.update = function(beta, sigmaSq, lambda2Sq) {
    p = length(beta)

    means = sqrt((lambda2Sq * sigmaSq) / (beta[2:p] - beta[1:p-1])^2)
    shape = lambda2Sq

    omegaSq.inv = numeric(p-1)
    for (i in 1:p-1) {
        omegaSq.inv[i] = rinvgauss(1, means[i], shape)
    }

    return(omegaSq.inv)
}

fused.sigmaSq.update = function(beta, tauSq.inv, omegaSq.inv, X, y) {
    p = length(beta)
    n = dim(X)[1]

    Sigma.beta.inv = calculate.Sigma.beta.inv(tauSq.inv, omegaSq.inv)

    v = y - X %*% beta
    rate = 0.5 * t(v) %*% v + 0.5 * t(beta) %*% Sigma.beta.inv %*% beta
    shape = (n - 1 + p) / 2

    sigmaSq.inv = rgamma(1, shape=shape, rate=rate)
    return(1 / sigmaSq.inv)
}

fused.lambda1Sq.update = function(p, r, tauSq.inv, delta) {
    tauSq = 1 / tauSq.inv
    shape = p + r
    rate = 0.5 * sum(tauSq^2) + delta
    lambda1Sq = rgamma(1, shape=shape, rate=rate)
    return(lambda1Sq)
}

fused.lambda2Sq.update = function(p, r, omegaSq.inv, delta) {
    omegaSq = 1 / omegaSq.inv
    shape = p + r - 1
    rate = 0.5 * sum(omegaSq^2) + delta
    lambda2Sq = rgamma(1, shape=shape, rate=rate)
    return(lambda2Sq)
}

fused.beta.estimate = function(X, y, r=1, delta=0.1, 
                        burnIn.iteration.count = 1000, 
                        sampling.iteration.count = 10000) {
    p = dim(X)[2]

    # Initializing parameters
    beta = rnorm(p)
    sigmaSq = rexp(1)
    tauSq.inv = rgamma(p, 1, 1)
    omegaSq.inv = rgamma(p-1, 1, 1)

    # Tuning parameters
    lambda1Sq = rgamma(1, 1, 0.1)
    lambda2Sq = rgamma(1, 1, 0.1)

    pb = txtProgressBar(0, burnIn.iteration.count + sampling.iteration.count, style=3)

    # Burn in Iterations
    for (iter in 1:burnIn.iteration.count) {
        beta = fused.beta.update(sigmaSq, tauSq.inv, omegaSq.inv, X, y)
        tauSq.inv = fused.tauSq.inv.update(beta, sigmaSq, lambda1Sq)
        omegaSq.inv = fused.omegaSq.inv.update(beta, sigmaSq, lambda2Sq)
        sigmaSq = fused.sigmaSq.update(beta, tauSq.inv, omegaSq.inv, X, y)
        lambda1Sq = fused.lambda1Sq.update(p, r, tauSq.inv, delta)
        lambda2Sq = fused.lambda2Sq.update(p, r, omegaSq.inv, delta)

        setTxtProgressBar(pb, iter)
    }

    beta.matrix = matrix(, nrow=sampling.iteration.count, ncol=p)
    sigmaSq.vector = numeric(length=sampling.iteration.count)

    # Sampling Iterations
    for (iter in 1:sampling.iteration.count) {
        beta = fused.beta.update(sigmaSq, tauSq.inv, omegaSq.inv, X, y)
        tauSq.inv = fused.tauSq.inv.update(beta, sigmaSq, lambda1Sq)
        omegaSq.inv = fused.omegaSq.inv.update(beta, sigmaSq, lambda2Sq)
        sigmaSq = fused.sigmaSq.update(beta, tauSq.inv, omegaSq.inv, X, y)
        lambda1Sq = fused.lambda1Sq.update(p, r, tauSq.inv, delta)
        lambda2Sq = fused.lambda2Sq.update(p, r, omegaSq.inv, delta)

        beta.matrix[iter,] = beta
        sigmaSq.vector[iter] = sigmaSq

        setTxtProgressBar(pb, iter + burnIn.iteration.count)
    }
    close(pb)

    return(colMeans(beta.matrix))
}

# N: number of simulations
gibbs.fused.lasso = function(X, true.beta, true.sigmaSq, N, r=1, delta=0.1, 
                            burnIn.iteration.count = 1000, sampling.iteration.count = 10000) {
    n = dim(X)[1]
    p = length(true.beta)
    beta.estimates = matrix(0, nrow=N, ncol=p)
    MSE.estimates = numeric(length=N)

    for (sim.i in 1:N) {
        print(paste0("Fused Lasso: ", N, " simulations, ", n, " observations"))
        print(paste0("Simulation ", sim.i, ":"))
        y = X %*% true.beta + rnorm(n, 0, sqrt(true.sigmaSq))

        beta = original.beta.estimate(X, y, r, delta, burnIn.iteration.count, sampling.iteration.count)
        MSE.estimates[sim.i] = t(y - X %*% beta) %*% (y - X %*% beta) / n
        beta.estimates[sim.i,] = beta

        clear_output()
    }

    beta.estimate = colMeans(beta.estimates)
    mse.estimate = mean(MSE.estimates)
    mse.std = sd(MSE.estimates)
    true.mse = colMeans((beta.estimates - matrix(rep(true.beta, N), ncol=p, byrow=TRUE))^2)

    return(list(
        beta.estimate=beta.estimate,
        mse.estimate=mse.estimate,
        true.beta=true.beta,
        true.mse=true.mse,
        mse.std=mse.std
    ))
}