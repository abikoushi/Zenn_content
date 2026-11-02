using Plots

# n次 Taylor 多項式
function poly_n_exp(x, n)
    return sum(x^k / factorial(k) for k in 0:n)
end

function myexp(x)
    log2 = 0.6931471805599453
    k, r = divrem(x, log2)
    return 2^k * poly_n_exp(r, 10)
end

plot(exp, 0, 4, label="exp")
plot!(myexp, label="myexp", linestyle=:dash)
png("myexp.png")

# GIF アニメーション
anim = @animate for n in 2:20
    x = range(0, 10, length=500)

    plot(
        x, exp.(x),
        label="exp(x)",
        linewidth=2,
        title="n = $n",
        xlabel="x",
        ylabel="y"
    )

    plot!(
        x,
        poly_n_exp.(x, n),
        label="Taylor polynomial",
        color=:orangered,
        linestyle=:dash,
        linewidth=2
    )
end

gif(anim, "taylor_exp.gif", fps=5)

#####

# n次 Taylor 多項式
function poly_n_log1p(x, n)
    return sum((-1)^(k+1) * x^k / k for k in 1:n)
end

function mylog(x)
    x > 0 || throw(DomainError(x))
    log2 = 0.6931471805599453
    m, e = frexp(x)
    r = 2m - 1.0
    k = e - 1
    return k*log2 + poly_n_log1p(r, 20)
end

plot(log, 0, 4, label="log")
plot!(mylog, label="mylog", linestyle=:dash)
png("mylog.png")

# GIF アニメーション
anim = @animate for n in 2:50
    x = range(0, 1.5, length=500)

    plot(
        x, log1p.(x),
        label="log(1+x)",
        linewidth=2,
        title="n = $n",
        xlabel="x",
        ylabel="y",
        ylims=(0.0, 1.0)
    )

    plot!(
        x,
        poly_n_log1p.(x, n),
        label="Taylor polynomial",
        color=:orangered,
        linestyle=:dash,
        linewidth=2
    )
end

gif(anim, "taylor_log1p.gif", fps=5)
