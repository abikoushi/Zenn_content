module vbgd

using Distributions
using Random
using Plots
using Zygote
using LinearAlgebra

function runGD(f, par0, lambda, lr, maxiter::Int) 
    rng = Random.default_rng()
    len = length(par0)
    mu = randn(rng,len)
    rho = randn(rng,len)
    logloss = zeros(maxiter)
    for i in 1:maxiter
        sigma = exp.(rho)
        epsilon = randn(rng,len)
        beta = mu + sigma.*epsilon
        gvec = f'(beta)
        g_mu = gvec + lambda*beta
        g_rho = (gvec.*epsilon + lambda*beta.*epsilon).*sigma .- 1
        mu -= lr * g_mu
        rho -= lr * g_rho
        logloss[i] = f(beta)
    end
    return mu, rho, logloss
end

function runNewton(f, par0, lambda, lr, maxiter::Int) 
    rng = Random.default_rng()
    len = length(par0)
    mu = randn(rng,len)
    rho = randn(rng,len)
    logloss = zeros(maxiter)
    for i in 1:maxiter
        sigma = exp.(rho)
        epsilon = randn(rng,len)
        beta = mu + sigma.*epsilon
        gvec = f'(beta)
        g2 = Zygote.hessian(beta0 -> f(beta0), beta)
        g_mu = gvec + lambda*beta
        g_rho = (gvec.*epsilon + lambda*beta.*epsilon).*sigma .- 1
        g2_mu = g2 + Diagonal(-lambda*sigma)
        g2_rho = g2 + Diagonal(sigma .*(epsilon .- lambda))
        mu -= lr  * inv(g2_mu) * g_mu
        rho -= lr * inv(g2_rho) * g_rho 
        logloss[i] = f(beta)
    end
    return mu, rho, logloss
end

end

using Distributions
using Random
using Plots
using Zygote
using LinearAlgebra

rng = Random.default_rng()
Random.seed!(1234)
x = sort(randn(rng,100))
X = [ones(100) x]
beta = [2.0, -1.0]

y = rand.(rng, Poisson.(exp.(X*beta)))

function poisonloss(beta,y,X)
    Xbeta = X*beta
    lambda = exp.(Xbeta)
    return -sum(y .* Xbeta - lambda)
end

f(beta) = poisonloss(beta,y,X)
    
betaini = [0.0,0.0]
@time μ, ρ, logloss = vbgd.runGD(f, betaini, 1.0, 1e-5, 100)
#計算時間は0.01秒くらい
plot(logloss,legend=false)

@time μ2, ρ2, logloss2 = vbgd.runNewton(f, betaini, 1.0, 0.001, 100)
#計算時間は0.01秒くらい
plot(logloss,legend=false)
plot!(logloss2,legend=false)

μ
μ2

eps = randn(2)
A = eps * eps'

inv(A) * eps
eps \ A

ε = randn(rng,2,1000)
betasmp = μ .+ exp.(ρ).*ε 

post = exp.(X*betasmp)
pred = rand.(rng,Poisson.(post))
predmean = mean(pred,dims=2)
lwr = [quantile(pred[i,:],0.025) for i in 1:100]
upr = [quantile(pred[i,:],0.975) for i in 1:100]

scatter(x,y,legend=false, tick_direction=:out, markerstrokewidth=0, color=:orange)
plot!(x, predmean, color="royalblue")
plot!(x,lwr,color="royalblue", linestyle = :dash)
plot!(x,upr,color="royalblue", linestyle = :dash)
#png("plot.png")