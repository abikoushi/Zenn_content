using Distributions
using Random
using Statistics
using Plots

N0 = Normal(0,1)
skewness(N0)
kurtosis(N0)

U1 = Uniform(-1,1)
skewness(U1)
kurtosis(U1)

L1 = Logistic(0,1)
skewness(L1)
kurtosis(L1)

p = plot(x -> pdf(N0,x), -4, 4, tick_direction=:out, label="Normal", linestyle=:solid)
plot!(p, x -> pdf(U1,x), label="Uniform", linestyle=:dash)
plot!(p, x -> pdf(L1,x), label="Logistic", linestyle=:dot)
png(p, "density1.png")

function simCLT(dist, size, iter)
    out = zeros(iter)
    rng = Random.default_rng()
    for i in 1:iter
        X = rand(rng, dist, size)
        out[i] = (mean(X) - mean(dist)) * sqrt(size)/std(dist)
    end
     return out
end

col1 = palette(:default)[1:3]
@time outU1 = simCLT(U1, 10, 10000)
h = histogram(outU1, legend=false, normalize=:pdf, color="lightgray", tick_direction=:out)
plot!(h, x->pdf(Normal(),x), linewidth=2, color=col1[1])
png(h, "hist1.png")


@time outL1 = simCLT(L1, 10, 10000)
h = histogram(outL1, legend=false, normalize=:pdf, color="lightgray", tick_direction=:out)
plot!(h, x->pdf(Normal(),x),  linewidth=2, color=col1[1])
png(h, "hist2.png")

E1 = Exponential(1)

skewness(E1) #2.0
kurtosis(E1) #6.0

p = plot(x -> pdf(N0,x), -4, 4, tick_direction=:out, label="Normal", linestyle=:solid)
plot!(x->pdf(E1,x), 0,4, label="Exponential", linestyle=:dash,tick_direction=:out)
png(p, "density2.png")



outE11 = simCLT(E1,  10, 10000)
h = histogram(outE11, legend=false, normalize=:pdf, color="lightgray")
plot!(h, x->pdf(Normal(),x), color="royalblue")
png(h, "hist3.png")

outE12 = simCLT(E1, 100, 10000)
h = histogram(outE12, legend=false, normalize=:pdf, color="lightgray")
plot!(h, x->pdf(Normal(),x), color="royalblue")
png(h, "hist4.png")
