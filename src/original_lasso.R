library(MASS)
library(statmod)
library(lars)
library(IRdisplay) # not necessary

original.beta.update = function(sigmaSq, tauSq.inv, X, y) {
    p = length(tauSq.inv)
    
    D.tauSq.inv = diag(tauSq.inv)
    
    A = solve(t(X) %*% X + D.tauSq.inv)

    mu0 = A %*% t(X) %*% y
    sigmaSq0 = sigmaSq * A

    beta = mvrnorm(1, mu=mu0, Sigma=sigmaSq0)
    return(beta)
}

original.tauSq.inv.update = function(beta, sigmaSq, lambdaSq) {
    p = length(beta)

    means = sqrt((lambdaSq * sqrt(sigmaSq)) / abs(beta))
    shape = lambdaSq

    tauSq.inv = numeric(p)
    for (i in 1:p) {
        tauSq.inv[i] = rinvgauss(1, means[i], shape)
    }

    return(tauSq.inv)
}


original.sigmaSq.update = function(beta, tauSq.inv, X, y) {
    p = length(beta)
    n = dim(X)[1]

    D.tauSq.inv = diag(tauSq.inv)

    v = y - X %*% beta
    rate = 0.5 * t(v) %*% v + 0.5 * t(beta) %*% D.tauSq.inv %*% beta
    shape = (n - 1 + p) / 2

    sigmaSq.inv = rgamma(1, shape=shape, rate=rate)
    return(1 / sigmaSq.inv)
}

original.lambdaSq.update = function(p, r, tauSq.inv, delta) {
    tauSq = 1 / tauSq.inv
    shape = p + r
    rate = 0.5 * sum(tauSq^2) + delta
    lambdaSq = rgamma(1, shape=shape, rate=rate)
    return(lambdaSq)
}

original.beta.estimate = function(X, y, r=1, delta=0.1, 
                        burnIn.iteration.count = 1000, 
                        sampling.iteration.count = 10000) {
    p = dim(X)[2]

    # Initializing parameters
    beta = rnorm(p)
    sigmaSq = rexp(1)
    tauSq.inv = rgamma(p, 1, 1)

    # Tuning parameters
    lambdaSq = rgamma(1, 1, 0.1)

    pb = txtProgressBar(0, burnIn.iteration.count + sampling.iteration.count, style=3)
 
    # Burn in Iterations
    for (iter in 1:burnIn.iteration.count) {
        beta = original.beta.update(sigmaSq, tauSq.inv, X, y)
        tauSq.inv = original.tauSq.inv.update(beta, sigmaSq, lambdaSq)
        sigmaSq = original.sigmaSq.update(beta, tauSq.inv, X, y)
        lambdaSq = original.lambdaSq.update(p, r, tauSq.inv, delta)

        setTxtProgressBar(pb, iter)
    }

    beta.matrix = matrix(, nrow=sampling.iteration.count, ncol=p)
    sigmaSq.vector = numeric(length=sampling.iteration.count)

    # Sampling Iterations
    for (iter in 1:sampling.iteration.count) {
        beta = original.beta.update(sigmaSq, tauSq.inv, X, y)
        tauSq.inv = original.tauSq.inv.update(beta, sigmaSq, lambdaSq)
        sigmaSq = original.sigmaSq.update(beta, tauSq.inv, X, y)
        lambdaSq = original.lambdaSq.update(p, r, tauSq.inv, delta)

        beta.matrix[iter,] = beta
        sigmaSq.vector[iter] = sigmaSq

        setTxtProgressBar(pb, iter + burnIn.iteration.count)
    }
    close(pb)

    return(colMeans(beta.matrix))
}

# N: Number of simulations
gibbs.original.lasso = function(X, true.beta, true.sigmaSq, N, r=1, delta=0.1, 
                                burnIn.iteration.count = 1000, sampling.iteration.count = 10000) {
    n = dim(X)[1]
    p = length(true.beta)
    beta.gibbs.estimates = matrix(0, nrow=N, ncol=p)
    MSE.gibbs.estimates = numeric(length=N)

    beta.lars.estimates = matrix(0, nrow=N, ncol=p)
    MSE.lars.estimates = numeric(length=N)

    for (sim.i in 1:N) {
        print(paste0("Original Lasso: ", N, " simulations, ", n, " observations"))
        print(paste0("Simulation ", sim.i, ":"))
        y = X %*% true.beta + rnorm(n, 0, sqrt(true.sigmaSq))

        # gibbs
        beta = original.beta.estimate(X, y, r, delta, burnIn.iteration.count, sampling.iteration.count)
        MSE.gibbs.estimates[sim.i] = t(y - X %*% beta) %*% (y - X %*% beta) / n
        beta.gibbs.estimates[sim.i,] = beta

        # lars
        beta = c(tail(coef(lars(X, y)), 1))
        MSE.lars.estimates[sim.i] = t(y - X %*% beta) %*% (y - X %*% beta) / n
        beta.lars.estimates[sim.i,] = beta

        clear_output() # not necessary
    }

    gibbs = list(
        beta.estimate = colMeans(beta.gibbs.estimates),
        mse.estimate = mean(MSE.gibbs.estimates),
        mse.std = sd(MSE.gibbs.estimates),
        true.mse = colMeans((beta.gibbs.estimates - matrix(rep(true.beta, N), ncol=p, byrow=TRUE))^2)
    )

    lars = list(
        beta.estimate = colMeans(beta.lars.estimates),
        mse.estimate = mean(MSE.lars.estimates),
        mse.std = sd(MSE.lars.estimates),
        true.mse = colMeans((beta.lars.estimates - matrix(rep(true.beta, N), ncol=p, byrow=TRUE))^2)
    )

    return(list(
        gibbs=gibbs, lars=lars, true.beta=true.beta
    ))
}