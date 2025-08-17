using Distributions
using Plots
using LogExpFunctions

p0 = repeat([1/6], 6)
-sum(xlogx, p0)

dice0 = Categorical(p0)
p = bar(x -> pdf(dice0, x), 1:6, ylims=[0,1], legend=false, color=:lightgrey, tick_direction=:out)
png(p, "bar_dice0.png")
entropy(dice0)

dice1 = Categorical( [0.1,0.05,0.1,0.2,0.15,0.4] )
p = bar(x -> pdf(dice1, x), 1:6, ylims=[0,1], legend=false, color=:lightgrey, tick_direction=:out)
entropy(dice1)
png(p, "bar_dice1.png")

dice2 = Categorical( [1,0,0,0,0,0] )
p = bar(x -> pdf(dice2, x), 1:6, ylims=[0,1.1], legend=false, color=:lightgrey, tick_direction=:out)
entropy(dice2)
png(p, "bar_dice2.png")


coin0 = Bernoulli(0.5)
p = plot(x -> kldivergence(coin0, Bernoulli(x)), legend=false)
xlabel!(p, "q")
ylabel!(p, "KLD")
png(p, "line_kld.png")

f(x, y) = kldivergence(Bernoulli(x), Bernoulli(y))
x = range(0.01, 0.99, length=100)
y = range(0.01, 0.99, length=100)
z = @. f(x', y)
p = heatmap(x, y, z)
xlabel!(p, "p")
ylabel!(p, "q")
png(p, "heat_kld.png")

