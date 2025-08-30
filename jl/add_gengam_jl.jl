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
export cdf, ccdf, pdf, logpdf, mean, loggeommean, shape, scale, power, params, rand
end

using StatsBase
using Statistics
using Plots
using Random
using RCall
using DataFrames

ggdist1 = GG.GeneralizedGamma(1,1,0.5)
ggdist2 = GG.GeneralizedGamma(1,2,2)
ggdist3 = GG.GeneralizedGamma(2,0.1,0.5)
ggdist4 = GG.GeneralizedGamma(2,1,2)
ggdist5 = GG.GeneralizedGamma(1,2.5,5)

R"
library(ggamma)
library(ggplot2)
library(patchwork)
dggamma2 = function(x,a,b,k){dggamma(x,b,k,a/k)}
"

x = 0.01:0.1:15
dens = GG.pdf(ggdist1, x)
df = DataFrame(x=x, y=dens)

R"
p1 = ggplot($df, aes(x, y))+
  geom_point()+
  stat_function(fun = dggamma2, args = list(a=1, b=1, k=0.5), n=1001)+
  ggtitle('(a=1, b=1, k=0.5)')+
  theme_bw()
"

dens = GG.pdf(ggdist2, x)
df = DataFrame(x=x, y=dens)

R"
p2 = ggplot($df, aes(x, y))+
  geom_point()+
  stat_function(fun = dggamma2, args = list(a=1, b=2, k=2))+
  ggtitle('(a=1, b=2, k=2)')+
  theme_bw()
"

dens = GG.pdf(ggdist3, x)
df = DataFrame(x=x, y=dens)

R"
p3 = ggplot($df, aes(x, y))+
  geom_point()+
  stat_function(fun = dggamma2, args = list(a=2, b=0.1, k=0.5))+
  ggtitle('(a=2, b=0.1, k=0.5)')+
  theme_bw()
"

dens = GG.pdf(ggdist4, x)
df = DataFrame(x=x, y=dens)

R"
p4 = ggplot($df, aes(x, y))+
  geom_point()+
  stat_function(fun = dggamma2, args = list(a=2, b=1, k=2))+
  ggtitle('(a=2, b=1, k=2)')+
  theme_bw()
"

dens = GG.pdf(ggdist5, x)
df = DataFrame(x=x, y=dens)

R"
p5 = ggplot($df, aes(x, y))+
  geom_point()+
  stat_function(fun = dggamma2, args = list(a=1, b=2.5, k=5))+
  ggtitle('(a=1, b=2.5, k=5)')+
  theme_bw()
"

R"
p = p1+p2+p3+p4+p5
ggsave('dggamma.png', plot=p, width=8, height=10)
"


##
rng1 = Random.default_rng(1)

x = GG.rand(rng1, ggdist1, 100000)
p1 = plot(y -> ecdf(x)(y), 0, 10, linestyle=:dash, label="ecdf")
plot!(p1, y -> GG.cdf(ggdist1,y), linestyle=:dashdot, label="true")
title!(p1, "(a=1,b=1,k=0.5)")

mean(x) #6.012540422658787
GG.mean(ggdist1) #6.0
mean(log, x) #0.848834299052846
GG.loggeommean(ggdist1) #0.8455686701969344

GG.rand!(rng1, ggdist2, x)
p2 = plot(y -> ecdf(x)(y), 0, 10, linestyle=:dash, label="ecdf")
plot!(p2, y -> GG.cdf(ggdist2,y), linestyle=:dashdot, label="true")
title!(p2, "(a=1,b=2,k=2)")

mean(x)
GG.mean(ggdist2)
mean(log, x)
GG.loggeommean(ggdist2)


GG.rand!(rng1, ggdist3, x)
p3 = plot(y -> ecdf(x)(y), 0, 10, linestyle=:dash, label="ecdf")
plot!(p3, y -> GG.cdf(ggdist3,y), linestyle=:dashdot, label="true")
title!(p3, "(a=2,b=0.1,k=0.5)")

mean(x) #1.9975355185928443
GG.mean(ggdist3) #2.0
mean(log, x)
GG.loggeommean(ggdist3)


GG.rand!(rng1, ggdist4, x)
p4 = plot(y -> ecdf(x)(y), 0, 10, linestyle=:dash, label="ecdf")
plot!(p4, y -> GG.cdf(ggdist4, y), linestyle=:dashdot, label="true")
title!(p4, "(a=2,b=1,k=2)")
mean(x)
GG.mean(ggdist4)
mean(log, x)
GG.loggeommean(ggdist4)

GG.rand!(rng1, ggdist5, x)
p5 = plot(y -> ecdf(x)(y), 0, 10, linestyle=:dash, label="ecdf")
plot!(p5, y -> GG.cdf(ggdist5, y), linestyle=:dashdot, label="true")
title!(p5, "(a=1,b=2.5,k=5)")

mean(x)
GG.mean(ggdist5)
mean(log, x)
GG.loggeommean(ggdist5)

p = plot(p1,p2,p3,p4,p5, size=(1000,800))
savefig("ecdf_ggamma.png")