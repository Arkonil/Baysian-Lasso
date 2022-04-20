library(MASS)
library(statmod)
library(IRdisplay)

group.beta.update = function(beta, groups, sigmaSq, tauSq.inv, X, y) {
    for (k in 1:length(groups)) {
        X.k = X[, groups[[k]]]
        A.k.inv = solve(t(X.k) %*% X.k + tauSq.inv[k]*diag(length(groups[[k]])))

        mu0 = A.k.inv %*% t(X[, groups[[k]]]) %*% (y - 0.5 * X[, -groups[[k]]] %*% beta[-groups[[k]]])
        sigmaSq0 = sigmaSq * A.k.inv

        beta[groups[[k]]] = mvrnorm(1, mu=mu0, Sigma=sigmaSq0)
    }
    return(beta)
}

group.tauSq.inv.update = function(beta, groups, sigmaSq, lambdaSq) {
    K = length(groups)

    tauSq.inv = numeric(K)
    for (k in 1:K) {
        mu0 = sqrt((lambdaSq * sigmaSq) / (t(beta[groups[[k]]]) %*% beta[groups[[k]]]))
        shape = lambdaSq

        tauSq.inv[k] = rinvgauss(1, mu0, shape)
    }

    return(tauSq.inv)
}

group.sigmaSq.update = function(beta, groups, tauSq.inv, X, y) {
    p = length(beta)
    n = dim(X)[1]

    v = y - X %*% beta
    rate = sum(v^2)
    for (k in 1:length(groups)) {
        rate = rate + tauSq.inv[k] * t(beta[groups[[k]]]) %*% beta[groups[[k]]]
    }
    rate = rate * 0.5    
    shape = (n - 1 + p) / 2

    sigmaSq.inv = rgamma(1, shape=shape, rate=rate)
    return(1 / sigmaSq.inv)
}

group.lambdaSq.update = function(p, r, tauSq.inv, delta) {
    K = length(tauSq.inv)
    tauSq = 1 / tauSq.inv
    shape = 0.5 * (p + K) + r
    rate = 0.5 * sum(tauSq^2) + delta
    lambdaSq = rgamma(1, shape=shape, rate=rate)
    return(lambdaSq)
}

group.beta.estimate = function(X, y, groups, r=1, delta=0.5, 
                            burnIn.iteration.count = 1000, sampling.iteration.count = 10000) {
    p = dim(X)[2]
    K = length(groups)

    # Initializing parameters
    beta = rnorm(p)
    sigmaSq = rexp(1)
    tauSq.inv = rgamma(K, 1, 1)

    # Tuning parameters
    lambdaSq = rgamma(1, 1, 0.1)

    pb = txtProgressBar(0, burnIn.iteration.count + sampling.iteration.count, style=3)

    # Burn in Iterations
    for (iter in 1:burnIn.iteration.count) {
        beta = group.beta.update(beta, groups, sigmaSq, tauSq.inv, X, y)
        tauSq.inv = group.tauSq.inv.update(beta, groups, sigmaSq, lambdaSq)
        sigmaSq = group.sigmaSq.update(beta, groups, tauSq.inv, X, y)
        lambdaSq = group.lambdaSq.update(p, r, tauSq.inv, delta)

        setTxtProgressBar(pb, iter)
    }

    beta.matrix = matrix(, nrow=sampling.iteration.count, ncol=p)
    sigmaSq.vector = numeric(length=sampling.iteration.count)

    # Sampling Iterations
    for (iter in 1:sampling.iteration.count) {
        beta = group.beta.update(beta, groups, sigmaSq, tauSq.inv, X, y)
        tauSq.inv = group.tauSq.inv.update(beta, groups, sigmaSq, lambdaSq)
        sigmaSq = group.sigmaSq.update(beta, groups, tauSq.inv, X, y)
        lambdaSq = group.lambdaSq.update(p, r, tauSq.inv, delta)

        beta.matrix[iter,] = beta
        sigmaSq.vector[iter] = sigmaSq

        setTxtProgressBar(pb, iter + burnIn.iteration.count)
    }
    close(pb)

    return(colMeans(beta.matrix))
}

gibbs.group.lasso = function(X, true.beta, true.sigmaSq, N, groups, r=1, delta=0.5, 
                            burnIn.iteration.count = 1000, sampling.iteration.count = 10000) {
    n = dim(X)[1]
    p = length(true.beta)

    beta.estimates = matrix(0, nrow=N, ncol=p)
    MSE.estimates = numeric(length=N)

    for (sim.i in 1:N) {
        print(paste0("Group Lasso: ", N, " simulations, ", n, " observations"))
        print(paste0("Simulation ", sim.i, ":"))
        y = X %*% true.beta + rnorm(n, 0, sqrt(true.sigmaSq))

        beta = group.beta.estimate(X, y, groups, r, delta, 
                            burnIn.iteration.count, sampling.iteration.count)
        MSE.estimates[sim.i] = t(y - X %*% beta) %*% (y - X %*% beta)
        beta.estimates[sim.i,] = beta

        clear_output()
    }

    beta.estimate = colMeans(beta.estimates)
    mse.estimate = sum(MSE.estimates) / (n * N)
    true.mse = colMeans((beta.estimates - matrix(rep(true.beta, N), ncol=p, byrow=TRUE))^2)

    return(list(
        beta.estimate=beta.estimate,
        mse.estimate=mse.estimate,
        true.beta=true.beta,
        true.mse=true.mse
    ))
}