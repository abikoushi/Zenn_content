module GG

using Random
using Distributions
using SpecialFunctions
import Distributions: ccdf, cdf, logpdf, pdf, quantile, mean, rand, params, shape, scale
import Distributions: @distr_support

struct GeneralizedGamma{Ta, Tb, Tk} <: ContinuousUnivariateDistribution
    a::Ta
    b::Tb
    k::Tk
end

@distr_support GeneralizedGamma 0.0 Inf

#### Parameters
params(d::GeneralizedGamma) = (d.a, d.b, d.k)
shape(d::GeneralizedGamma) = d.a
scale(d::GeneralizedGamma) = d.b
power(d::GeneralizedGamma) = d.k

#### Evaluation
function gamma_cdf(a, x)
    return gamma_inc(a,x,0)[1]
end

function gamma_ccdf(a, x)
    return gamma_inc(a,x,0)[2]
end


function cdf(d::GeneralizedGamma, x::Real)
    if x < zero(x)
        return zero(x)
    else
        shp, scl, pwr = params(d)
        return gamma_cdf(shp/pwr, (x/scl)^pwr) 
    end
end

function ccdf(d::GeneralizedGamma, x::Real)
    if x < zero(x)
        return one(x)
    else
        shp, scl, pwr = params(d)
        return gamma_ccdf(shp/pwr, (x/scl)^pwr)
    end
end

function pdf(d::GeneralizedGamma, x::Real)
    if x < zero(x)
        return zero(x)
    else
        shp, scl, pwr = params(d)
        return (pwr/(scl^shp))*x^(shp-1)*exp(-(x/scl)^pwr)/gamma(shp/pwr)
    end
end

function logpdf(d::GeneralizedGamma, x::Real)
    if x < zero(x)
        return zero(x)
    else
        shp, scl, pwr = params(d)
        return log(pwr)-shp*log(scl) + (shp-1)*log(x) - (x/scl)^pwr -loggamma(shp/pwr)
    end
end

function quantile(d::GeneralizedGamma, p)
    shp, scl, pwr = params(d)
    r = gamma_inc_inv(shp/pwr, p, 1-p)
    return scl*(r^inv(pwr))
end

function rand(rng::AbstractRNG, d::GeneralizedGamma)
    p = rand(rng)
    return quantile(d, p)
end


#### Statistics

function mean(d::GeneralizedGamma)
    shp, scl, pwr = params(d)
    return scl*gamma((shp+1)/pwr)/gamma(shp/pwr)
end

function loggeommean(d::GeneralizedGamma)
    shp, scl, pwr = params(d)
    return digamma(shp/pwr)/pwr + log(scl)
end

export GeneralizedGamma
export cdf, ccdf, pdf, logpdf, mean, meanlog, shape, scale, power, params, rand
end

using Statistics
using Plots
using Random
using RCall

R"
library(ggamma)
curve(dggamma(x,1,1,0.5),0,8)
"

plot(x-> GG.pdf(GG.GeneralizedGamma(1,1,0.5), x),0, 8, label = "(1,1,0.5)")
plot!(x-> GG.pdf(GG.GeneralizedGamma(1,2,2), x), label = "(1,2,2)")
plot!(x-> GG.pdf(GG.GeneralizedGamma(2,1,0.5), x), label = "(2,1,0.5)")
plot!(x-> GG.pdf(GG.GeneralizedGamma(2,1,2), x), label = "(2,1,2)")
plot!(x-> GG.pdf(GG.GeneralizedGamma(1,5,5), x), label = "(1,5,5)")
png("gr.png")

rng1 = Random.default_rng(1)
x = GG.rand(rng1,GG.GeneralizedGamma(1,1,0.5),100000)
mean(x)
GG.mean(GG.GeneralizedGamma(1,1,0.5))
mean(log, x)
GG.loggeommean(GG.GeneralizedGamma(1,1,0.5))

GG.rand!(rng1,GG.GeneralizedGamma(1,2,2), x)
mean(x)
GG.mean(GG.GeneralizedGamma(1,2,2))
mean(log, x)
GG.loggeommean(GG.GeneralizedGamma(1,2,2))

GG.rand!(MersenneTwister(1),GG.GeneralizedGamma(2,1,0.5), x)
mean(x)
GG.mean(GG.GeneralizedGamma(2,1,0.5))
mean(log, x)
GG.meanlog(GG.GeneralizedGamma(2,1,0.5))

GG.rand!(MersenneTwister(1),GG.GeneralizedGamma(2,1,2), x)
mean(x)
GG.mean(GG.GeneralizedGamma(2,1,2))
mean(log, x)
GG.meanlog(GG.GeneralizedGamma(2,1,2))

GG.rand!(MersenneTwister(1),GG.GeneralizedGamma(1,5,5), x)
mean(x)
GG.mean(GG.GeneralizedGamma(1,5,5))
mean(log, x)
GG.meanlog(GG.GeneralizedGamma(1,5,5))
