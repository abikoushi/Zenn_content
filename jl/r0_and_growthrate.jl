using DifferentialEquations
using DataFrames
using CairoMakie, AlgebraOfGraphics
using Distributions

function SIR!(du, u, p, t)
    S, I, R = u
    β, γ = p
    du[1] = -β * S * I
    du[2] = β * S * I - γ * I
    du[3] = γ * I
end

function solve_SIR(β, γ, S0, I0, R0, tspan)
    u0 = [S0, I0, R0]
    p = [β, γ]
    prob = ODEProblem(SIR!, u0, tspan, p)
    sol = solve(prob, dt = 0.1, adaptive = false)
    return sol
end


β = 1
γ = 0.25
S0 = 0.999
I0 = 0.001
R0 = 0.0
tspan = (0.0, 100.0)
solution = solve_SIR(β, γ, S0, I0, R0,  tspan)

solution_df = DataFrame(time=solution.t, S=solution[1, :], I=solution[2, :], R=solution[3, :])

p1 = AlgebraOfGraphics.data(solution_df) * 
    mapping(:time, :I) *
    visual(Lines)
p2 = AlgebraOfGraphics.data(solution_df)*
    mapping(:time, :S) *
    visual(Lines, color=:lightgray)
p3 = AlgebraOfGraphics.data(solution_df)*
    mapping(:time, :R) *
    visual(Lines, color=:lightgray)

p3 + p2 + p1 |> draw(axis=(ylabel="I",))

exp_df = DataFrame(time=solution.t,  exp = I0*exp.((β-γ) .* solution.t))

p4 = AlgebraOfGraphics.data(exp_df)*
    mapping(:time, :exp) *
    visual(Lines, color=:blue, linestyle=:dash) 

fig = draw(p4+p3+p2+p1,
 axis = (ylabel="I", limits = (nothing, (-0.1,1.1))),
 figure=(size=(600, 500),) 
)

save("plot1.png", fig)

#####

function R0hat_gamma(m, CV)
    r = 0.135
    ab2 = (m * CV)^2
    b = ab2 / m
    a = m / b
    val = mgf(Gamma(a, b), -r)
    return 1 / val
end

means = 1:0.1:6
CVs   = 0.1:0.1:1.0

resdf = DataFrame(
    mean = Float64[],
    CV   = Float64[],
    R0   = Float64[]
)

for m in means, cv in CVs
    push!(resdf, (
        mean = m,
        CV   = cv,
        R0   = R0hat_gamma(m, cv)
    ))
end


plt_heat = data(resdf) *
    mapping(:mean, :CV, :R0) *
    visual(Heatmap)

fig2 = draw(plt_heat,  figure=(size=(500, 500),))

save("plot2.png", fig2)