using Distributions
using Random
using Statistics
using Plots

G1 = Gamma(7,1)
X = rand(G1, 100)
Y = X .+ 10
p = histogram(X, legend=false, alpha=0.5)
histogram!(p, Y, alpha=0.5)
plot(p)

G2 = Gamma(0.8,1)
skewness(G1)
kurtosis(G1)
skewness(G2)
kurtosis(G2)

function simCLT(dist, iter, size)
    out = zeros(iter)
    rng = Random.default_rng()
    for i in 1:iter
        X = rand(rng, dist, size)
        out[i] = mean(X .- mean(dist)) * sqrt(size)/std(dist)
    end
     return out
end

out1 = simCLT(G1, 10000, 10)
h = histogram(out1, legend=false, normalize=:pdf, color="lightgray")
plot!(h, x->pdf(Normal(),x), color="royalblue")

out2 = simCLT(G2, 10000, 10)
h = histogram(out2, legend=false, normalize=:pdf, color="lightgray")
plot!(h, x->pdf(Normal(),x), color="royalblue")
