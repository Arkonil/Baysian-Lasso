library(MASS)
library(statmod)
library(elasticnet)
library(IRdisplay)

calculate.D.tauSq.inv = function(tauSq.inv, lambda2) {
    return(diag(1 / (tauSq.inv + lambda2)))
}

elasticnet.beta.update = function(sigmaSq, tauSq.inv, lambda2, X, y) {
    p = length(tauSq.inv)
    
    D.tauSq.inv = calculate.D.tauSq.inv(tauSq.inv, lambda2)
    
    A = solve(t(X) %*% X + D.tauSq.inv)

    mu0 = A %*% t(X) %*% y
    sigmaSq0 = sigmaSq * A

    beta = mvrnorm(1, mu=mu0, Sigma=sigmaSq0)
    return(beta)
}

elasticnet.tauSq.inv.update = function(beta, sigmaSq, lambda1Sq) {
    p = length(beta)

    means = sqrt((lambda1Sq * sigmaSq) / beta^2)
    shape = lambda1Sq

    tauSq.inv = numeric(p)
    for (i in 1:p) {
        tauSq.inv[i] = rinvgauss(1, means[i], shape)
    }

    return(tauSq.inv)
}

elasticnet.sigmaSq.update = function(beta, tauSq.inv, lambda2, X, y) {
    p = length(beta)
    n = dim(X)[1]

    D.tauSq.inv = calculate.D.tauSq.inv(tauSq.inv, lambda2)

    v = y - X %*% beta
    rate = 0.5 * t(v) %*% v + 0.5 * t(beta) %*% D.tauSq.inv %*% beta
    shape = (n - 1 + p) / 2

    sigmaSq.inv = rgamma(1, shape=shape, rate=rate)
    return(1 / sigmaSq.inv)
}

elasticnet.lambda1Sq.update = function(p, r1, tauSq.inv, delta1) {
    tauSq = 1 / tauSq.inv
    shape = p + r1
    rate = 0.5 * sum(tauSq^2) + delta1
    lambda1Sq = rgamma(1, shape=shape, rate=rate)
    return(lambda1Sq)
}

elasticnet.lambda2.update = function(p, r2, beta, sigmaSq, delta2) {
    shape = 0.5 * p + r2
    rate = sum(beta^2) / (2 * sigmaSq) + delta2
    lambda2 = rgamma(1, shape=shape, rate=rate)
    return(lambda2)
}

elasticnet.beta.estimate = function(X, y, r1=1, r2=1, delta1=0.1, delta2=0.1, 
                        burnIn.iteration.count = 1000, 
                        sampling.iteration.count = 10000) {
    p = dim(X)[2]

    # Initializing parameters
    beta = rnorm(p)
    sigmaSq = rexp(1)
    tauSq.inv = rgamma(p, 1, 1)

    # Tuning parameters
    lambda1Sq = rgamma(1, 1, 0.1)
    lambda2 = rgamma(1, 1, 0.1)

    pb = txtProgressBar(0, burnIn.iteration.count + sampling.iteration.count, style=3)

    # Burn in Iterations
    for (iter in 1:burnIn.iteration.count) {
        beta = elasticnet.beta.update(sigmaSq, tauSq.inv, lambda2, X, y)
        tauSq.inv = elasticnet.tauSq.inv.update(beta, sigmaSq, lambda1Sq)
        sigmaSq = elasticnet.sigmaSq.update(beta, tauSq.inv, lambda2, X, y)
        lambda1Sq = elasticnet.lambda1Sq.update(p, r1, tauSq.inv, delta1)
        lambda2 = elasticnet.lambda2.update(p, r2, beta, sigmaSq, delta2)

        setTxtProgressBar(pb, iter)
    }

    beta.matrix = matrix(, nrow=sampling.iteration.count, ncol=p)
    sigmaSq.vector = numeric(length=sampling.iteration.count)

    # Sampling Iterations
    for (iter in 1:sampling.iteration.count) {
        beta = elasticnet.beta.update(sigmaSq, tauSq.inv, lambda2, X, y)
        tauSq.inv = elasticnet.tauSq.inv.update(beta, sigmaSq, lambda1Sq)
        sigmaSq = elasticnet.sigmaSq.update(beta, tauSq.inv, lambda2, X, y)
        lambda1Sq = elasticnet.lambda1Sq.update(p, r1, tauSq.inv, delta1)
        lambda2 = elasticnet.lambda2.update(p, r2, beta, sigmaSq, delta2)

        beta.matrix[iter,] = beta
        sigmaSq.vector[iter] = sigmaSq

        setTxtProgressBar(pb, iter + burnIn.iteration.count)
    }
    close(pb)

    return(colMeans(beta.matrix))
}

# N: number of simulations
gibbs.elasticnet.lasso = function(X, true.beta, true.sigmaSq, N, r1=1, r2=1, delta1=0.1, delta2=0.1, 
                                burnIn.iteration.count = 1000, sampling.iteration.count = 10000) {
    n = dim(X)[1]
    p = length(true.beta)
    beta.gibbs.estimates = matrix(0, nrow=N, ncol=p)
    MSE.gibbs.estimates = numeric(length=N)

    beta.lars.estimates = matrix(0, nrow=N, ncol=p)
    MSE.lars.estimates = numeric(length=N)

    for (sim.i in 1:N) {
        print(paste0("ElasticNet Lasso: ", N, " simulations, ", n, " observations"))
        print(paste0("Simulation ", sim.i, ":"))
        y = X %*% true.beta + rnorm(n, 0, sqrt(true.sigmaSq))

        # gibbs
        beta = elasticnet.beta.estimate(X, y, r1, r2, delta1, delta2, burnIn.iteration.count, sampling.iteration.count)
        MSE.gibbs.estimates[sim.i] = t(y - X %*% beta) %*% (y - X %*% beta) / n
        beta.gibbs.estimates[sim.i,] = beta

        # lars
        beta = c(tail(predict(enet(X, y))$coefficients, 1))
        MSE.lars.estimates[sim.i] = t(y - X %*% beta) %*% (y - X %*% beta) / n
        beta.lars.estimates[sim.i,] = beta

        clear_output()
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