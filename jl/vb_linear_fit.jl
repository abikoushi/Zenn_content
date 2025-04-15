module vblm
function vb_linear_fit(X, y, a0=1e-2, b0=1e-4, c0=1e-2, d0=1e-4)
    # Prior parameters
    N, D = size(X)
    X_corr = X' * X
    Xy_corr = X' * y
    an = a0 + N / 2
    gammaln_an = lgamma(an)
    cn = c0 + D / 2
    gammaln_cn = lgamma(cn)
    # Iterate to find hyperparameters
    L_last = -Inf
    max_iter = 500
    E_a = c0 / d0
    for iter in 1:max_iter
        # Covariance and weight of linear model
        invV = E_a * I(D) + X_corr
        chol_V = cholesky(invV)
        V = inv(chol_V) * inv(chol_V)'
        logdetV = -2 * sum(log(diag(chol_V.U)))
        w = V * Xy_corr
        # Parameters of noise model (an remains constant)
        sse = sum((X * w - y) .^ 2)
        bn = b0 + 0.5 * (sse + E_a * (w' * w))
        E_t = an / bn
        # Hyperparameters of covariance prior (cn remains constant)
        dn = d0 + 0.5 * (E_t * (w' * w) + tr(V))
        E_a = cn / dn
        # Variational bound, ignoring constant terms for now
        L = -0.5 * (E_t * sse + sum(X .* (X * V))) + 0.5 * logdetV - b0 * E_t \
            + gammaln_an - an * log(bn) + an + gammaln_cn - cn * log(dn)
        # Variational bound must grow
        if L_last > L
            println("Last bound ", L_last, ", current bound ", L)
            error("Variational bound should not reduce")
        end
        # Stop if change in variation bound is < 0.001%
        if abs(L_last - L) < abs(0.00001 * L)
            break
        end
        L_last = L
    end
    if iter == max_iter
        @warn "Bayesian linear regression reached maximum number of iterations."
    end
    # Augment variational bound with constant terms
    L -= 0.5 * (N * log(2 * pi) - D) - lgamma(a0) + a0 * log(b0) - lgamma(c0) + c0 * log(d0)
    return w, V, invV, logdetV, an, bn, E_a, L
end
end

